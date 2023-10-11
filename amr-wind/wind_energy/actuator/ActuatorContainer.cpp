#include "amr-wind/wind_energy/actuator/ActuatorContainer.H"
#include "amr-wind/wind_energy/actuator/Actuator.H"
#include "amr-wind/core/gpu_utils.H"
#include "amr-wind/core/Field.H"

#include "AMReX_Scan.H"

#include <AMReX_Print.H>
#include <algorithm>

namespace amr_wind::actuator {

ActuatorCloud::ActuatorCloud(const int nobjects)
    : num_pts(nobjects, 0), global_id(nobjects, -1), num_objects(nobjects)
{}

ActuatorContainer::ActuatorContainer(
    amrex::AmrCore& mesh, const int num_objects)
    : amrex::AmrParticleContainer<
          NumPStructReal,
          NumPStructInt,
          NumPArrayReal,
          NumPArrayInt>(&mesh)
    , m_mesh(mesh)
    , m_data(num_objects)
    , m_proc_pos(amrex::ParallelDescriptor::NProcs(), vs::Vector::zero())
    , m_pos_device(amrex::ParallelDescriptor::NProcs(), vs::Vector::zero())
    , m_proc_offsets(amrex::ParallelDescriptor::NProcs() + 1, 0)
    , m_proc_offsets_device(amrex::ParallelDescriptor::NProcs() + 1)
{}

/** Allocate memory and initialize the particles within the container
 *
 *  This method is only called once during the simulation. It allocates the
 *  arrays for holding the position vector and velocity data on host memory and
 *  also initializes corresponding particles within the first available particle
 *  tile within the container. It is expected that the actuator manager instance
 *  has already populated the number of points per turbine before invoking this
 *  method.
 */
void ActuatorContainer::initialize_container()
{
    BL_PROFILE("amr-wind::actuator::ActuatorContainer::initialize_container");

    compute_local_coordinates();

    // Initialize global data arrays
    const int total_pts =
        std::accumulate(m_data.num_pts.begin(), m_data.num_pts.end(), 0);
    m_data.position.resize(total_pts);
    m_data.velocity.resize(total_pts);
    m_data.density.resize(total_pts);
    m_data.est_force.resize(total_pts);

    {
        const int nproc = amrex::ParallelDescriptor::NProcs();
        amrex::Vector<int> pts_per_proc(nproc, 0);
#ifdef AMREX_USE_MPI
        MPI_Allgather(
            &total_pts, 1, MPI_INT, pts_per_proc.data(), 1, MPI_INT,
            amrex::ParallelDescriptor::Communicator());
#else
        pts_per_proc[0] = total_pts;
#endif
        m_proc_offsets[0] = 0;
        for (int i = 1; i <= nproc; ++i) {
            m_proc_offsets[i] = m_proc_offsets[i - 1] + pts_per_proc[i - 1];
        }
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_proc_offsets.begin(),
            m_proc_offsets.end(), m_proc_offsets_device.begin());
    }

    initialize_particles(total_pts);
}

void ActuatorContainer::initialize_particles(const int total_pts)
{
    // Initialize particle container data structures.
    //
    // We assign all particles into the first available container within this
    // MPI rank and let redistribute take care of scattering the particles into
    // the respective MPI rank containing the cell enclosing this particle. The
    // actual redistribution happens after position vectors are updated in
    // update_positions.

    // query the particle ID from the container. We should always be starting
    // from 1.
    ParticleType::NextID(1U);
    const auto id_start = ParticleType::NextID();
    AMREX_ALWAYS_ASSERT(id_start == 1U);
    const int iproc = amrex::ParallelDescriptor::MyProc();

    // Flag indicating if a tile was found where all particles were deposited.
    bool assigned = false;
    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; ((lev < nlevels) && !assigned); ++lev) {
        for (auto mfi = MakeMFIter(lev); (mfi.isValid() && !assigned); ++mfi) {
            auto& ptile = GetParticles(
                lev)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
            AMREX_ASSERT(ptile.empty());
            ptile.resize(total_pts);
            auto* pstruct = ptile.GetArrayOfStructs()().data();

            amrex::ParallelFor(
                total_pts, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
                    auto& pp = pstruct[ip];

                    pp.id() = id_start + ip;
                    pp.cpu() = iproc;
                    pp.idata(0) = ip;
                });
            assigned = true;
        }
    }

    // Indicate that we have initialized the containers and remaining methods
    // are safe to use
    m_container_initialized = true;
    m_is_scattered = false;
}

void ActuatorContainer::reset_container()
{
    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; lev < nlevels; ++lev) {
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto* pstruct = pti.GetArrayOfStructs()().data();

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
                auto& pp = pstruct[ip];
                pp.id() = -1;
            });
        }
    }
    Redistribute();

    const int total_pts = static_cast<int>(m_data.velocity.size());
    initialize_particles(total_pts);
}

/** Update position vectors of the particles within a container based on the
 *  data provided by actuator instances.
 *
 *  This method assumes that the particles have been restored to their
 *  originating MPI ranks after any scattering. This happens during
 *  initialization (after AcutatorContainer::initialize_particles) and after
 *  ActuatorContainer::sample_velocities.
 *
 *  After executing this method, the particles have been scattered across the
 *  domain such that each particle is contained within a tile that has the cell
 *  enclosing the particle.
 *
 */
void ActuatorContainer::update_positions()
{
    BL_PROFILE("amr-wind::actuator::ActuatorContainer::update_positions");
    AMREX_ALWAYS_ASSERT(m_container_initialized && !m_is_scattered);

    const auto dpos = gpu::device_view(m_data.position);
    const auto* const dptr = dpos.data();
    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; lev < nlevels; ++lev) {
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto* pstruct = pti.GetArrayOfStructs()().data();

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
                auto& pp = pstruct[ip];
                const auto idx = pp.idata(0);

                const auto& pvec = dptr[idx];
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    pp.pos(n) = pvec[n];
                }
            });
        }
    }

    // Scatter particles to appropriate MPI ranks
    Redistribute();

    // Indicate that it is safe to sample velocities
    m_is_scattered = true;
}

/** Interpolate the velocity field using a trilinear interpolation
 *
 *  This method performs three tasks:
 *    - Sample the velocity field and interpolate onto particle location
 *    - Restore the particles back to their original MPI rank
 *    - Copy data from particles into the buffer used by actuator instances
 *
 *  After invocation if this method, the particles are back in their original
 *  rank and the position vectors can be updated safely.
 */
void ActuatorContainer::sample_fields(const Field& vel, const Field& density)
{
    BL_PROFILE("amr-wind::actuator::ActuatorContainer::sample_velocities");
    AMREX_ALWAYS_ASSERT(m_container_initialized && m_is_scattered);

    // Sample velocity field
    interpolate_fields(vel, density);

    // Recall particles to the MPI ranks that contains their corresponding
    // turbines
    // Redistribute();

    // Populate the velocity buffer that all actuator instances can access
    populate_field_buffers();

    // Indicate that the particles have been restored to their original MPI rank
    m_is_scattered = false;
}

/** Helper method for ActuatorContainer::sample_fields
 *
 *  Loops over the particle tiles and copies velocity data from particles to the
 *  velocity and density arrays.
 */
void ActuatorContainer::populate_field_buffers()
{
    BL_PROFILE("amr-wind::actuator::ActuatorContainer::populate_vel_buffer");
    const size_t num_buff_entries =
        m_proc_offsets.back() * static_cast<size_t>(NumPStructReal);

    amrex::Vector<amrex::Real> buff_host(num_buff_entries);
    amrex::Gpu::DeviceVector<amrex::Real> buff_device(num_buff_entries, 0.0);

    auto* buffer_pointer = buff_device.data();
    auto* offsets = m_proc_offsets_device.data();

    const int nlevels = m_mesh.finestLevel() + 1;

    for (int lev = 0; lev < nlevels; ++lev) {
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto* pstruct = pti.GetArrayOfStructs()().data();

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
                auto& pp = pstruct[ip];
                const auto iproc = pp.cpu();
                const auto idx = offsets[iproc] + pp.idata(0);

                for (int n = 0; n < NumPStructReal; ++n) {
                    buffer_pointer[idx * NumPStructReal + n] = pp.rdata(n);
                }
            });
        }
    }

    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, buff_device.begin(), buff_device.end(),
        buff_host.begin());
#ifdef AMREX_USE_MPI
    const int num_entires = static_cast<int>(buff_host.size());
    MPI_Allreduce(
        MPI_IN_PLACE, buff_host.data(), num_entires, MPI_DOUBLE, MPI_SUM,
        amrex::ParallelDescriptor::Communicator());
#endif
    {
        auto& vel_arr = m_data.velocity;
        auto& den_arr = m_data.density;
        auto& est_arr = m_data.est_force;
        const int npts = static_cast<int>(vel_arr.size());
        const int ioff = m_proc_offsets[amrex::ParallelDescriptor::MyProc()];
        for (int i = 0; i < npts; ++i) {
            // See NumPStructReal def in ActuatorContainer.H
            for (int j = 0; j < AMREX_SPACEDIM; ++j) {
                vel_arr[i][j] = buff_host[(ioff + i) * NumPStructReal + j];
            }
            den_arr[i] =
                buff_host[(ioff + i) * NumPStructReal + AMREX_SPACEDIM];
            for (int j = 0; j < AMREX_SPACEDIM; ++j) {
                est_arr[i][j] = buff_host[(ioff + i) * NumPStructReal + AMREX_SPACEDIM+1 + j];
            }
        }
    }
}

/** Helper method for ActuatorContainer::sample_fields
 *
 *  Performs a trilinear interpolation of the velocity/desnity field to particle
 *  locations. It also updates the particle locations such that the next
 *  Redistribute call restores the particles back to their original MPI rank
 *  where they were created.
 */
void ActuatorContainer::interpolate_fields(
    const Field& vel, const Field& density)
{
    BL_PROFILE("amr-wind::actuator::ActuatorContainer::interpolate_velocities");
    auto* dptr = m_pos_device.data();
    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = m_mesh.Geom(lev);
        const auto dx = geom.CellSizeArray();
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();

        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto* pstruct = pti.GetArrayOfStructs()().data();
            const auto varr = vel(lev).const_array(pti);
            const auto darr = density(lev).const_array(pti);

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
                auto& pp = pstruct[ip];
                // Determine offsets within the containing cell
                const amrex::Real x =
                    (pp.pos(0) - plo[0] - 0.5 * dx[0]) * dxi[0];
                const amrex::Real y =
                    (pp.pos(1) - plo[1] - 0.5 * dx[1]) * dxi[1];
                const amrex::Real z =
                    (pp.pos(2) - plo[2] - 0.5 * dx[2]) * dxi[2];

                // Index of the low corner
                const int i = static_cast<int>(std::floor(x));
                const int j = static_cast<int>(std::floor(y));
                const int k = static_cast<int>(std::floor(z));

                // Interpolation weights in each direction (linear basis)
                const amrex::Real wx_hi = (x - i);
                const amrex::Real wy_hi = (y - j);
                const amrex::Real wz_hi = (z - k);

                const amrex::Real wx_lo = 1.0 - wx_hi;
                const amrex::Real wy_lo = 1.0 - wy_hi;
                const amrex::Real wz_lo = 1.0 - wz_hi;

                const int iproc = pp.cpu();

                //amrex::AllPrint() << "[iproc="<<iproc<<"]"
                //    << " particle " << pp.idata(0)
                //    << " ("<<pp.pos(0)<<","<<pp.pos(1)<<","<<pp.pos(2)<<")"
                //    << " --> ("<<x<<","<<y<<","<<z<<")"
                //    << " at ("<<i<<","<<j<<","<<k<<")"
                //    << std::endl;

                // velocity
                for (int ic = 0; ic < AMREX_SPACEDIM; ++ic) {
                    pp.rdata(ic) =
                        wx_lo * wy_lo * wz_lo * varr(i, j, k, ic) +
                        wx_lo * wy_lo * wz_hi * varr(i, j, k + 1, ic) +
                        wx_lo * wy_hi * wz_lo * varr(i, j + 1, k, ic) +
                        wx_lo * wy_hi * wz_hi * varr(i, j + 1, k + 1, ic) +
                        wx_hi * wy_lo * wz_lo * varr(i + 1, j, k, ic) +
                        wx_hi * wy_lo * wz_hi * varr(i + 1, j, k + 1, ic) +
                        wx_hi * wy_hi * wz_lo * varr(i + 1, j + 1, k, ic) +
                        wx_hi * wy_hi * wz_hi * varr(i + 1, j + 1, k + 1, ic);

                    // Reset position vectors so that the particles return back
                    // to the MPI ranks with the turbines upon redistribution
                    pp.pos(ic) = dptr[iproc][ic];
                }

                // density
                pp.rdata(AMREX_SPACEDIM) =
                    wx_lo * wy_lo * wz_lo * darr(i, j, k) +
                    wx_lo * wy_lo * wz_hi * darr(i, j, k + 1) +
                    wx_lo * wy_hi * wz_lo * darr(i, j + 1, k) +
                    wx_lo * wy_hi * wz_hi * darr(i, j + 1, k + 1) +
                    wx_hi * wy_lo * wz_lo * darr(i + 1, j, k) +
                    wx_hi * wy_lo * wz_hi * darr(i + 1, j, k + 1) +
                    wx_hi * wy_hi * wz_lo * darr(i + 1, j + 1, k) +
                    wx_hi * wy_hi * wz_hi * darr(i + 1, j + 1, k + 1);
            });
        }
    }
}

/** Identify and mark actuator surface faces
 *
 *  Used by the ThinBody actuator. Based on interpolate_fields().
 */
void ActuatorContainer::mark_surface_faces(
    IntField& xface,
    IntField& yface,
    IntField& zface,
    IntField& mask_cell,
    IntField& mask_node)
{
    BL_PROFILE("amr-wind::actuator::ActuatorContainer::mark_surface_faces");
    AMREX_ALWAYS_ASSERT(m_container_initialized && m_is_scattered);

    const bool DEBUG = false;

    if (DEBUG)
        amrex::Print() << "Marking surface faces" << std::endl;
    xface.setVal(0);
    yface.setVal(0);
    zface.setVal(0);
    mask_cell.setVal(0);
    mask_node.setVal(0);

    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = m_mesh.Geom(lev);
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();

        auto& xffab = xface(lev);
        auto& yffab = yface(lev);
        auto& zffab = zface(lev);
        auto& mcfab = mask_cell(lev);
        auto& mnfab = mask_node(lev);

        // Note: This particle iterator code will only execute on grids
        // _that contain particles_
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto* pstruct = pti.GetArrayOfStructs()().data();
            const auto& xf_arr = xffab.array(pti);
            const auto& yf_arr = yffab.array(pti);
            const auto& zf_arr = zffab.array(pti);
            const auto& mc_arr = mcfab.array(pti);
            const auto& mn_arr = mnfab.array(pti);

            if (DEBUG)
                amrex::Print() << "pti (lev="<<lev<<")"
                    << " Array4("<<yf_arr.begin.x<<":"<<yf_arr.end.x-1<<","
                                 <<yf_arr.begin.y<<":"<<yf_arr.end.y-1<<","
                                 <<yf_arr.begin.z<<":"<<yf_arr.end.z-1<<")"
                    << std::endl;

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
                auto& pp = pstruct[ip];
                // Determine offsets _from the nearest node_
                const amrex::Real x = (pp.pos(0) - plo[0]) * dxi[0];
                const amrex::Real y = (pp.pos(1) - plo[1]) * dxi[1];
                const amrex::Real z = (pp.pos(2) - plo[2]) * dxi[2];

                // Index of the low corner
                constexpr amrex::Real domain_eps = 1.0e-6;
                const int i = static_cast<int>(std::floor(x + domain_eps));
                const int j = static_cast<int>(std::floor(y + domain_eps));
                const int k = static_cast<int>(std::floor(z + domain_eps));

                if (DEBUG)
                    amrex::AllPrint() << "[iproc="<<pp.cpu()<<"]"
                        << " mark particle " << pp.idata(0)
                        << " ("<<pp.pos(0)<<","<<pp.pos(1)<<","<<pp.pos(2)<<")"
                        << " --> ("<<x<<","<<y<<","<<z<<")"
                        << " at ("<<i<<","<<j<<","<<k<<")"
                        << std::endl;

                // x,y,z will be integers if the requested position is
                // staggered, unlike interpolate_fields. If a surface face is
                // found, mark the face and the adjacent cells/nodes. Each face
                // will have 2 adjacent cells and 4 adjacent nodes.
                if (x==i) {
                    xf_arr(i  ,j  ,k  ) = -1; // staggered in i
                    mc_arr(i  ,j  ,k  ) = -1;
                    mc_arr(i-1,j  ,k  ) = -1;
                    mn_arr(i  ,j  ,k  ) = -1;
                    mn_arr(i  ,j+1,k  ) = -1;
                    mn_arr(i  ,j  ,k+1) = -1;
                    mn_arr(i  ,j+1,k+1) = -1;
                }
                if (y==j) {
                    yf_arr(i  ,j  ,k  ) = -1; // staggered in j
                    mc_arr(i  ,j  ,k  ) = -1;
                    mc_arr(i  ,j-1,k  ) = -1;
                    mn_arr(i  ,j  ,k  ) = -1;
                    mn_arr(i+1,j  ,k  ) = -1;
                    mn_arr(i  ,j  ,k+1) = -1;
                    mn_arr(i+1,j  ,k+1) = -1;
                }
                if (z==k) {
                    zf_arr(i  ,j  ,k  ) = -1; // staggered in k
                    mc_arr(i  ,j  ,k  ) = -1;
                    mc_arr(i  ,j  ,k-1) = -1;
                    mn_arr(i  ,j  ,k  ) = -1;
                    mn_arr(i+1,j  ,k  ) = -1;
                    mn_arr(i  ,j+1,k  ) = -1;
                    mn_arr(i+1,j+1,k  ) = -1;
                }
            });
        }

        // At this point, faces/cells/nodes will be set to -1, but we still
        // need to handle the fringe case in which a face/cell lies on grid
        // boundaries or ghost cells/nodes. The particle iteration will only
        // select one of the adjacent grids (assuming the actuator point lies
        // on a face) and therefore not all grids will necessarily have the
        // same face/cells/nodes flagged.
        xffab.SumBoundary(geom.periodicity());
        yffab.SumBoundary(geom.periodicity());
        zffab.SumBoundary(geom.periodicity());
        mcfab.SumBoundary(geom.periodicity());
        mnfab.SumBoundary(geom.periodicity());

        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto* pstruct = pti.GetArrayOfStructs()().data();
            const auto& mc_arr = mcfab.array(pti);
            const auto& mn_arr = mnfab.array(pti);

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
                auto& pp = pstruct[ip];
                // Determine offsets _from the nearest node_
                const amrex::Real x = (pp.pos(0) - plo[0]) * dxi[0];
                const amrex::Real y = (pp.pos(1) - plo[1]) * dxi[1];
                const amrex::Real z = (pp.pos(2) - plo[2]) * dxi[2];

                // Index of the low corner
                constexpr amrex::Real domain_eps = 1.0e-6;
                const int i = static_cast<int>(std::floor(x + domain_eps));
                const int j = static_cast<int>(std::floor(y + domain_eps));
                const int k = static_cast<int>(std::floor(z + domain_eps));

                // If grids overlap at more than one edge or face, SumBoundary
                // will end up with mask values < -1. Need to make sure fringes
                // are exactly -1 at this point, otherwise the resulting mask
                // will have unexpected or undefined values.
                if (mc_arr(i  ,j  ,k  ) < -1) mc_arr(i  ,j  ,k  ) = -1;
                if (mc_arr(i-1,j  ,k  ) < -1) mc_arr(i-1,j  ,k  ) = -1;
                if (mc_arr(i  ,j-1,k  ) < -1) mc_arr(i  ,j-1,k  ) = -1;
                if (mc_arr(i  ,j  ,k-1) < -1) mc_arr(i  ,j  ,k-1) = -1;
                if (mn_arr(i  ,j  ,k  ) < -1) mn_arr(i  ,j  ,k  ) = -1;
                if (mn_arr(i  ,j+1,k  ) < -1) mn_arr(i  ,j+1,k  ) = -1;
                if (mn_arr(i  ,j  ,k+1) < -1) mn_arr(i  ,j  ,k+1) = -1;
                if (mn_arr(i  ,j+1,k+1) < -1) mn_arr(i  ,j+1,k+1) = -1;
                if (mn_arr(i+1,j  ,k  ) < -1) mn_arr(i+1,j  ,k  ) = -1;
                if (mn_arr(i  ,j  ,k+1) < -1) mn_arr(i  ,j  ,k+1) = -1;
                if (mn_arr(i+1,j  ,k+1) < -1) mn_arr(i+1,j  ,k+1) = -1;
                if (mn_arr(i+1,j  ,k  ) < -1) mn_arr(i+1,j  ,k  ) = -1;
                if (mn_arr(i  ,j+1,k  ) < -1) mn_arr(i  ,j+1,k  ) = -1;
                if (mn_arr(i+1,j+1,k  ) < -1) mn_arr(i+1,j+1,k  ) = -1;

                if ((mc_arr(i  ,j  ,k  ) > 0) ||
                    (mc_arr(i-1,j  ,k  ) > 0) ||
                    (mc_arr(i  ,j-1,k  ) > 0) ||
                    (mc_arr(i  ,j  ,k-1) > 0) ||
                    (mn_arr(i  ,j  ,k  ) > 0) ||
                    (mn_arr(i  ,j+1,k  ) > 0) ||
                    (mn_arr(i  ,j  ,k+1) > 0) ||
                    (mn_arr(i  ,j+1,k+1) > 0) ||
                    (mn_arr(i+1,j  ,k  ) > 0) ||
                    (mn_arr(i  ,j  ,k+1) > 0) ||
                    (mn_arr(i+1,j  ,k+1) > 0) ||
                    (mn_arr(i+1,j  ,k  ) > 0) ||
                    (mn_arr(i  ,j+1,k  ) > 0) ||
                    (mn_arr(i+1,j+1,k  ) > 0))
                    amrex::Abort("Mask should be either -1 or 0 at this point");
            });
        }

        // Now all surface faces (and adjacent cells/nodes) will be -1, so add
        // 1 to have the masks be:
        //   0 on (and adjacent to) surface faces
        //   1 everywhere else
        int nghost = mask_cell.num_grow()[0];
        xffab.plus(1,0,1); // val, comp, num_comp, nghost=0
        yffab.plus(1,0,1);
        zffab.plus(1,0,1);
        mcfab.plus(1,0,1,nghost);
        mnfab.plus(1,0,1,nghost);
    }
}

/** Estimate force based on a local control volume and no penetration condition
 *
 *  Used by the ThinBody actuator. Based on mark_surface_faces(). Momentum sum
 *  occurs during the ICNS advective update.
 */
void ActuatorContainer::estimate_force(
    const Field& mom_sum)
{
    BL_PROFILE("amr-wind::actuator::ActuatorContainer::estimate_force");
    AMREX_ALWAYS_ASSERT(m_container_initialized && m_is_scattered);

    const bool DEBUG = false;

    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = m_mesh.Geom(lev);
        const auto dx = geom.CellSizeArray();
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();

        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto* pstruct = pti.GetArrayOfStructs()().data();
            const auto cv_arr = mom_sum(lev).array(pti);

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
                auto& pp = pstruct[ip];
                // Determine offsets _from the nearest node_
                const amrex::Real x = (pp.pos(0) - plo[0]) * dxi[0];
                const amrex::Real y = (pp.pos(1) - plo[1]) * dxi[1];
                const amrex::Real z = (pp.pos(2) - plo[2]) * dxi[2];

                // Index of the low corner
                constexpr amrex::Real domain_eps = 1.0e-6;
                const int i = static_cast<int>(std::floor(x + domain_eps));
                const int j = static_cast<int>(std::floor(y + domain_eps));
                const int k = static_cast<int>(std::floor(z + domain_eps));

                // x,y,z will be integers if the requested position is
                // staggered, unlike interpolate_fields
                //
                // See NumPStructReal def in ActuatorContainer.H
                const int idx = AMREX_SPACEDIM + 1;
                amrex::Real xcc = plo[0] + (i+0.5)*dx[0];
                amrex::Real ycc = plo[1] + (j+0.5)*dx[1];
                amrex::Real zcc = plo[2] + (k+0.5)*dx[2];
                if (x==i) {
                    // surface is on a x-face
                    pp.rdata(idx  ) = cv_arr(i-1,j,k,0) + cv_arr(i,j,k,0);
                    pp.rdata(idx+1) = cv_arr(i-1,j,k,1) + cv_arr(i,j,k,1);
                    pp.rdata(idx+2) = cv_arr(i-1,j,k,2) + cv_arr(i,j,k,2);
                    if (DEBUG) {
                        amrex::AllPrint() << "[iproc="<<pp.cpu()<<"]"
                            << " particle " << pp.idata(0)
                            << " ("<<pp.pos(0)<<","<<pp.pos(1)<<","<<pp.pos(2)<<")"
                            << " xsurf i,j,k = "<<i<<","<<j<<","<<k
                            << " est_force = ("<<cv_arr(i-1,j,k,0)<<","<<cv_arr(i-1,j,k,1)<<","<<cv_arr(i-1,j,k,2)<<")"
                            <<             " ("<<cv_arr(i  ,j,k,0)<<","<<cv_arr(i  ,j,k,1)<<","<<cv_arr(i  ,j,k,2)<<")"
                            //<< " at ("<<xcc-dx[0]<<","<<ycc<<","<<zcc<<")"
                            //<<    " ("<<xcc      <<","<<ycc<<","<<zcc<<")"
                            << std::endl;
                        AMREX_ALWAYS_ASSERT(std::abs((xcc-dx[0]) - (pp.pos(0)-dx[0]/2)) < 1e-14);
                        AMREX_ALWAYS_ASSERT(std::abs( xcc        - (pp.pos(0)+dx[0]/2)) < 1e-14);
                        AMREX_ALWAYS_ASSERT(std::abs( ycc        -  pp.pos(1)         ) < 1e-14);
                        AMREX_ALWAYS_ASSERT(std::abs( zcc        -  pp.pos(2)         ) < 1e-14);
                    }
                }
                if (y==j) {
                    // surface is on a y-face
                    pp.rdata(idx  ) = cv_arr(i,j-1,k,0) + cv_arr(i,j,k,0);
                    pp.rdata(idx+1) = cv_arr(i,j-1,k,1) + cv_arr(i,j,k,1);
                    pp.rdata(idx+2) = cv_arr(i,j-1,k,2) + cv_arr(i,j,k,2);
                    if (DEBUG) {
                        amrex::AllPrint() << "[iproc="<<pp.cpu()<<"]"
                            << " particle " << pp.idata(0)
                            << " ("<<pp.pos(0)<<","<<pp.pos(1)<<","<<pp.pos(2)<<")"
                            << " ysurf i,j,k = "<<i<<","<<j<<","<<k
                            << " est_force = ("<<cv_arr(i,j-1,k,0)<<","<<cv_arr(i,j-1,k,1)<<","<<cv_arr(i,j-1,k,2)<<")"
                            <<             " ("<<cv_arr(i,j  ,k,0)<<","<<cv_arr(i,j  ,k,1)<<","<<cv_arr(i,j  ,k,2)<<")"
                            //<< " at ("<<xcc<<","<<ycc-dx[1]<<","<<zcc<<")"
                            //<<    " ("<<xcc<<","<<ycc      <<","<<zcc<<")"
                            << std::endl;
                        AMREX_ALWAYS_ASSERT(std::abs( xcc        -  pp.pos(0)         ) < 1e-14);
                        AMREX_ALWAYS_ASSERT(std::abs((ycc-dx[1]) - (pp.pos(1)-dx[1]/2)) < 1e-14);
                        AMREX_ALWAYS_ASSERT(std::abs( ycc        - (pp.pos(1)+dx[1]/2)) < 1e-14);
                        AMREX_ALWAYS_ASSERT(std::abs( zcc        -  pp.pos(2)         ) < 1e-14);
                    }
                }
                if (z==k) {
                    // surface is on a z-face
                    pp.rdata(idx  ) = cv_arr(i,j,k-1,0) + cv_arr(i,j,k,0);
                    pp.rdata(idx+1) = cv_arr(i,j,k-1,1) + cv_arr(i,j,k,1);
                    pp.rdata(idx+2) = cv_arr(i,j,k-1,2) + cv_arr(i,j,k,2);
                    if (DEBUG) {
                        amrex::AllPrint() << "[iproc="<<pp.cpu()<<"]"
                            << " particle " << pp.idata(0)
                            << " ("<<pp.pos(0)<<","<<pp.pos(1)<<","<<pp.pos(2)<<")"
                            << " zsurf i,j,k = "<<i<<","<<j<<","<<k
                            << " est_force = ("<<cv_arr(i,j,k-1,0)<<","<<cv_arr(i,j,k-1,1)<<","<<cv_arr(i,j,k-1,2)<<")"
                            <<             " ("<<cv_arr(i,j,k  ,0)<<","<<cv_arr(i,j,k  ,1)<<","<<cv_arr(i,j,k  ,2)<<")"
                            //<< " at ("<<xcc<<","<<ycc<<","<<zcc-dx[2]<<")"
                            //<<    " ("<<xcc<<","<<ycc<<","<<zcc      <<")"
                            << std::endl;
                        AMREX_ALWAYS_ASSERT(std::abs( xcc        -  pp.pos(0)         ) < 1e-14);
                        AMREX_ALWAYS_ASSERT(std::abs( ycc        -  pp.pos(1)         ) < 1e-14);
                        AMREX_ALWAYS_ASSERT(std::abs((zcc-dx[2]) - (pp.pos(2)-dx[2]/2)) < 1e-14);
                        AMREX_ALWAYS_ASSERT(std::abs( zcc        - (pp.pos(2)+dx[2]/2)) < 1e-14);
                    }
                }
            });
        }
    }
}

/** Determine position vector of a point within each MPI rank
 *
 *  Loops over the boxArray and determines the first patch that belongs to the
 *  current MPI rank. Uses that to generate a position vector that is known to
 *  exist in a given rank. Setting the position vectors of the particles to this
 *  know location will recall particles belonging to this rank during
 *  Redistribute.
 */
void ActuatorContainer::compute_local_coordinates()
{
    BL_PROFILE(
        "amr-wind::actuator::ActuatorContainer::compute_local_coordinates");
    const int nprocs = amrex::ParallelDescriptor::NProcs();
    const int iproc = amrex::ParallelDescriptor::MyProc();

    // Reset position vectors to zero (required for parallel reduce sum)
    m_proc_pos.assign(nprocs, vs::Vector::zero());

    // Flag indicating whether a point within the domain belonging to this rank
    // has been found
    bool assigned = false;
    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; ((lev < nlevels) && !assigned); ++lev) {
        const auto& ba = m_mesh.boxArray(lev);
        const auto& dm = m_mesh.DistributionMap(lev);

        const int nbx = static_cast<int>(dm.size());
        for (int i = 0; (i < nbx) && !assigned; ++i) {
            if (dm[i] != iproc) {
                continue;
            }

            const auto& geom = m_mesh.Geom(lev);
            const auto& bx = ba[i];
            const int* lo = bx.loVect();

            auto& pvec = m_proc_pos[iproc];
            pvec.x() = geom.ProbLo()[0] + (lo[0] + 0.5) * geom.CellSize()[0];
            pvec.y() = geom.ProbLo()[1] + (lo[1] + 0.5) * geom.CellSize()[1];
            pvec.z() = geom.ProbLo()[2] + (lo[2] + 0.5) * geom.CellSize()[2];

            // Indicate that we have found a point and it is safe to exit the
            // loop
            assigned = true;
        }
    }

    // Share position vectors with every process
    amrex::ParallelDescriptor::ReduceRealSum(
        &(m_proc_pos[0].x()),
        static_cast<int>(m_proc_pos.size()) * vs::Vector::ncomp);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_proc_pos.begin(), m_proc_pos.end(),
        m_pos_device.begin());
}

} // namespace amr_wind::actuator
