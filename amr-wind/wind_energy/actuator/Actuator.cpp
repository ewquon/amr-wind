#include "amr-wind/wind_energy/actuator/Actuator.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "amr-wind/wind_energy/actuator/ActuatorContainer.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldRepo.H"

#include <algorithm>
#include <memory>

namespace amr_wind::actuator {

Actuator::Actuator(CFDSim& sim)
    : m_sim(sim), m_act_source(sim.repo().declare_field("actuator_src_term", 3))
{}

Actuator::~Actuator() = default;

void Actuator::pre_init_actions()
{
    BL_PROFILE("amr-wind::actuator::Actuator::post_init_actions");
    amrex::ParmParse pp(identifier());

    amrex::Vector<std::string> labels;
    pp.getarr("labels", labels);

    pp.query("weighted_sampling", m_weighted_sampling);

    pp.query("limit_source_CFL", m_limit_source_cfl);
    pp.query("CFL_limit", m_cfl_limit);

    const int nturbines = static_cast<int>(labels.size());

    for (int i = 0; i < nturbines; ++i) {
        const std::string& tname = labels[i];
        const std::string& prefix = identifier() + "." + tname;
        amrex::ParmParse pp1(prefix);

        std::string type;
        pp.query("type", type);
        pp1.query("type", type);
        AMREX_ALWAYS_ASSERT(!type.empty());

        auto obj = ActuatorModel::create(type, m_sim, tname, i);

        const std::string default_prefix = identifier() + "." + type;
        utils::ActParser inp(default_prefix, prefix);

        obj->read_inputs(inp);
        m_actuators.emplace_back(std::move(obj));
    }
}

void Actuator::post_init_actions()
{
    BL_PROFILE("amr-wind::actuator::Actuator::post_init_actions");

    amrex::Vector<int> act_proc_count(amrex::ParallelDescriptor::NProcs(), 0);
    for (auto& act : m_actuators) {
        act->determine_root_proc(act_proc_count);
    }

    {
        // Sanity check that we have processed the turbines correctly
        int nact =
            std::accumulate(act_proc_count.begin(), act_proc_count.end(), 0);
        AMREX_ALWAYS_ASSERT(num_actuators() == nact);
    }

    for (auto& act : m_actuators) {
        act->init_actuator_source();
    }

    setup_container();
    update_positions();
    update_velocities();
    compute_forces();
    compute_source_term();
    prepare_outputs();
}

void Actuator::post_regrid_actions()
{
    for (auto& act : m_actuators) {
        act->determine_influenced_procs();
    }

    setup_container();
}

void Actuator::pre_advance_work()
{
    BL_PROFILE("amr-wind::actuator::Actuator::pre_advance_work");

    m_container->reset_container();
    update_positions();
    update_velocities();
    compute_forces();
    compute_source_term();
    communicate_turbine_io();
}

void Actuator::communicate_turbine_io()
{
#ifdef AMR_WIND_USE_HELICS
    if (!m_sim.helics().is_activated()) {
        return;
    }
    // send power and yaw from root actuator proc to io proc
    const int ptag = 0;
    const int ytag = 1;
    const size_t size = 1;
    for (auto& ac : m_actuators) {
        if (ac->info().is_root_proc) {
            amrex::ParallelDescriptor::Send(
                &m_sim.helics().m_turbine_power_to_controller[ac->info().id],
                size, amrex::ParallelDescriptor::IOProcessorNumber(), ptag);
            amrex::ParallelDescriptor::Send(
                &m_sim.helics()
                     .m_turbine_wind_direction_to_controller[ac->info().id],
                size, amrex::ParallelDescriptor::IOProcessorNumber(), ytag);
        }
        if (amrex::ParallelDescriptor::IOProcessor()) {
            amrex::ParallelDescriptor::Recv(
                &m_sim.helics().m_turbine_power_to_controller[ac->info().id],
                size, ac->info().root_proc, ptag);
            amrex::ParallelDescriptor::Recv(
                &m_sim.helics()
                     .m_turbine_wind_direction_to_controller[ac->info().id],
                size, ac->info().root_proc, ytag);
        }
    }
#endif
}

/** Set up the container for sampling velocities
 *
 *  Allocates memory and initializes the particles corresponding to actuator
 *  nodes for all turbines that influence the current MPI rank. This method is
 *  invoked once during initialization and during regrid step.
 */
void Actuator::setup_container()
{
    const int ntotal = num_actuators();
    const int nlocal = static_cast<int>(std::count_if(
        m_actuators.begin(), m_actuators.end(),
        [](const std::unique_ptr<ActuatorModel>& obj) {
            return obj->info().sample_vel_in_proc;
        }));

    m_container = std::make_unique<ActuatorContainer>(m_sim.mesh(), nlocal);

    auto& pinfo = m_container->m_data;
    for (int i = 0, il = 0; i < ntotal; ++i) {
        if (m_actuators[i]->info().sample_vel_in_proc) {
            pinfo.global_id[il] = i;
            pinfo.num_pts[il] = m_actuators[i]->num_velocity_points(); // == num_span_pts
            ++il;
        }
    }

    m_container->initialize_container();
}

/** Update actuator positions and sample velocities at new locations.
 *
 *  This method loops over all the turbines local to this MPI rank and updates
 *  the position vectors. These new locations are provided to the sampling
 *  container that samples velocities at these new locations.
 *
 *  \sa Actuator::update_velocities
 */
void Actuator::update_positions()
{
    BL_PROFILE("amr-wind::actuator::Actuator::update_positions");
    auto& pinfo = m_container->m_data;
    for (int i = 0, ic = 0; i < pinfo.num_objects; ++i) {
        const auto ig = pinfo.global_id[i];
        auto vpos =
            ::amr_wind::utils::slice(pinfo.position, ic, pinfo.num_pts[i]);
        m_actuators[ig]->update_positions(vpos);
        ic += pinfo.num_pts[i];
    }

    if (!m_weighted_sampling) // -- original approach w/ particles
    {
        m_container->update_positions();

        // Mark ThinBody actuator surface faces
        if (m_sim.repo().field_exists("mom_flux_sum")) {
            auto& xface = m_sim.repo().get_int_field("mom_xface_mask");
            auto& yface = m_sim.repo().get_int_field("mom_yface_mask");
            auto& zface = m_sim.repo().get_int_field("mom_zface_mask");
            m_container->mark_surface_faces(xface,yface,zface);
        }

        // Sample velocities at the new locations
        const auto& vel = m_sim.repo().get_field("velocity");
        const auto& density = m_sim.repo().get_field("density");
        m_container->sample_fields(vel, density);
    }
    else
    {
        amrex::Print() << "[Actuator::update_positions] Sampling weighted velocities (DEPRECATED)" << std::endl;
        sample_weighted_velocities();
    }
}

/** Helper method to calculate the average velocities at actuator grid
 *  points weighted by the body-force projection function
 */
void Actuator::sample_weighted_velocities()
{
    BL_PROFILE("amr-wind::actuator::Actuator::sample_weighted_velocities");
    const int nlevels = m_sim.repo().num_active_levels();
    const auto& vfield = m_sim.repo().get_field("velocity");

    //for (int lev = 0; lev < nlevels; ++lev) {
    // assume the finest level captures most of the projection function
    int lev = nlevels - 1;
    {
        const auto& geom = m_sim.mesh().Geom(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(vfield(lev)); mfi.isValid(); ++mfi) {
            for (auto& ac : m_actuators) {
                //if (ac->info().sample_vel_in_proc)
                ac->integrate_weighted_velocities(lev, mfi, geom);
            }
        }
    }
}

/** Provide updated velocities from container to actuator instances
 *
 *  \sa Acuator::update_positions
 */
void Actuator::update_velocities()
{
    BL_PROFILE("amr-wind::actuator::Actuator::update_velocities");

    auto& pinfo = m_container->m_data;
    for (int i = 0, ic = 0; i < pinfo.num_objects; ++i) {
        const auto ig = pinfo.global_id[i];

        const auto vel =
            ::amr_wind::utils::slice(pinfo.velocity, ic, pinfo.num_pts[i]);

        const auto density =
            ::amr_wind::utils::slice(pinfo.density, ic, pinfo.num_pts[i]);

        if (m_weighted_sampling) {
            // This will just call ops::UpdateVelOp
            m_actuators[ig]->update_fields(vel, density, false);
        } else {
            // This will copy the particle velocities to the actuator grid
            // prior to calling ops::UpdateVelOp
            m_actuators[ig]->update_fields(vel, density, true);
        }
        ic += pinfo.num_pts[i];
    }
}

/** Helper method to compute forces on all actuator components
 */
void Actuator::compute_forces()
{
    BL_PROFILE("amr-wind::actuator::Actuator::compute_forces");
    for (auto& ac : m_actuators) {
        if (ac->info().actuator_in_proc) {
            ac->compute_forces();
        }
    }
}

void Actuator::compute_source_term()
{
    BL_PROFILE("amr-wind::actuator::Actuator::compute_source_term");
    m_act_source.setVal(0.0);
    const int nlevels = m_sim.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        auto& sfab = m_act_source(lev);
        const auto& geom = m_sim.mesh().Geom(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(sfab); mfi.isValid(); ++mfi) {
            for (auto& ac : m_actuators) {
                if (ac->info().actuator_in_proc) {
                    ac->compute_source_term(lev, mfi, geom);
                }
            }
        }
    }

    if (m_limit_source_cfl)
    {
        const int finest_level = m_sim.mesh().finestLevel();
        auto& sfab_fine = m_act_source(finest_level);
        auto const& src_arr = sfab_fine.const_arrays();
        amrex::Real srcmax_lev = 0.0;
        srcmax_lev += amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpMax>{},
            amrex::TypeList<amrex::Real>{},
            sfab_fine,
            amrex::IntVect(0),
            [=] AMREX_GPU_HOST_DEVICE(
                int box_no, int i, int j, int k) -> amrex::GpuTuple<amrex::Real> {
                auto const& s_bx = src_arr[box_no];
                return amrex::max<amrex::Real>(
                    amrex::Math::abs(s_bx(i, j, k, 0)),
                    amrex::Math::abs(s_bx(i, j, k, 1)),
                    amrex::Math::abs(s_bx(i, j, k, 2)),
                    static_cast<amrex::Real>(-1.0));
            });
        amrex::ParallelAllReduce::Max<amrex::Real>(
                    srcmax_lev, amrex::ParallelContext::CommunicatorSub());
        //amrex::Print(amrex::Print::AllProcs) << "max src val = " << srcmax_lev << std::endl;

        // source CFL = (src_term / (rho*dx))**0.5 * dt
        // max_source_mag = rho * dx * CFL_limit**2 / dt**2  [force/volume]
        const auto& dx = m_sim.mesh().Geom(finest_level).CellSize();
        const amrex::Real dt = m_sim.time().deltaT();
        const amrex::Real rho = 1.0;
        auto max_source_mag_x = rho * dx[0] * m_cfl_limit * m_cfl_limit / (dt*dt);
        auto max_source_mag_y = rho * dx[1] * m_cfl_limit * m_cfl_limit / (dt*dt);
        auto max_source_mag_z = rho * dx[2] * m_cfl_limit * m_cfl_limit / (dt*dt);
        amrex::Real gain = amrex::min<amrex::Real>(
            max_source_mag_x / srcmax_lev,
            max_source_mag_y / srcmax_lev,
            max_source_mag_z / srcmax_lev,
            1.0 // don't scale forces if current srcmax_lev < max_source_mag
        );
        amrex::Print() << "max source magnitude = " << srcmax_lev << std::endl;
        amrex::Print() << "max allowable source magnitude = "
            << max_source_mag_x << " "
            << max_source_mag_y << " "
            << max_source_mag_z << std::endl;
        amrex::Print() << "gain = " << gain << std::endl;

        // scale source term field
        for (int lev = 0; lev < nlevels; ++lev) {
            auto& sfab = m_act_source(lev);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(sfab); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                const auto& sarr = sfab.array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    sarr(i, j, k, 0) *= gain;
                    sarr(i, j, k, 1) *= gain;
                    sarr(i, j, k, 2) *= gain;
                });
            }
        }

        // scale actuator forces (point forces, integrated lift/drag)
        for (auto& ac : m_actuators) {
            if (ac->info().actuator_in_proc) {
                ac->scale_forces(gain);
            }
        }
    }
}

void Actuator::prepare_outputs()
{
    const std::string out_dir_prefix = "post_processing/actuator";
    const std::string sname =
        amrex::Concatenate(out_dir_prefix, m_sim.time().time_index());
    if (!amrex::UtilCreateDirectory(sname, 0755)) {
        amrex::CreateDirectoryFailed(sname);
    }
    const int iproc = amrex::ParallelDescriptor::MyProc();
    for (auto& ac : m_actuators) {
        if (ac->info().root_proc == iproc) {
            ac->prepare_outputs(sname);
        }
    }
}

void Actuator::post_advance_work()
{
    const int iproc = amrex::ParallelDescriptor::MyProc();
    for (auto& ac : m_actuators) {
        if (ac->info().root_proc == iproc) {
            ac->write_outputs();
        }
    }
}

} // namespace amr_wind::actuator
