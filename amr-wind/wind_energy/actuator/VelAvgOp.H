#ifndef VELAVGOP_H_
#define VELAVGOP_H_

#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/wind_energy/actuator/actuator_utils.H"
#include "amr-wind/core/FieldRepo.H"

namespace amr_wind::actuator::ops {

template <typename ActTrait>
class VelAvgOp
{
private:
    typename ActTrait::DataType& m_data;
    const Field& m_vfield;
  //const Field& m_rfield;

    DeviceVecList m_pos;
    DeviceVecList m_vel;
    DeviceVecList m_epsilon;
    DeviceTensorList m_orientation;

    amrex::Gpu::DeviceVector<amrex::Real> m_sum;
    amrex::Vector<amrex::Real> m_host_sum;

    void copy_to_device();
    void copy_to_host(); // update ActGrid.vel

public:
    explicit VelAvgOp(typename ActTrait::DataType& data)
        : m_data(data)
        , m_vfield(m_data.sim().repo().get_field("velocity"))
      //, m_rfield(m_data.sim().repo().get_field("density"))
    {}

    void initialize();

    void setup_op() {
        copy_to_device();

        // reset velocities before integrating
        const int npts = m_pos.size();
        for (int ip = 0; ip < npts; ++ip) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                m_sum[ip] = 0.0;
                m_vel[ip][idim] = 0.0;
            }
        }
    }

    void post();

    // integrate weighted velocities -- on device
    void operator()(
        const int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom);
};

template <typename ActTrait>
void VelAvgOp<ActTrait>::initialize()
{
    const auto& grid = m_data.grid();
    m_pos.resize(grid.pos.size());
    m_vel.resize(grid.vel.size());
    m_epsilon.resize(grid.epsilon.size());
    m_orientation.resize(grid.orientation.size());
    m_sum.resize(grid.pos.size());
    m_host_sum.resize(grid.pos.size());
}

template <typename ActTrait>
void VelAvgOp<ActTrait>::copy_to_device()
{
    const auto& grid = m_data.grid();

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.pos.begin(), grid.pos.end(),
        m_pos.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.epsilon.begin(), grid.epsilon.end(),
        m_epsilon.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.orientation.begin(),
        grid.orientation.end(), m_orientation.begin());
}

template <typename ActTrait>
void VelAvgOp<ActTrait>::copy_to_host()
{
    auto& grid = m_data.grid();

    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, m_vel.begin(), m_vel.end(),
        grid.vel.begin());

    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, m_sum.begin(), m_sum.end(),
        m_host_sum.begin());
}

template <typename ActTrait>
void VelAvgOp<ActTrait>::operator()(
    const int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom)
{
    const std::string fname = ActTrait::identifier();
    BL_PROFILE("amr-wind::VelAvgOp<" + fname + ">");

    const auto& bx = mfi.tilebox();
    const auto& varr = m_vfield(lev).array(mfi);
  //const auto& rarr = m_rfield(lev).array(mfi);
    const auto& problo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();

    amrex::Real vol = dx[0];
    for (int idim = 1; idim < AMREX_SPACEDIM; ++idim)
        vol *= dx[idim];

    const int npts = m_pos.size();
    const auto* pos = m_pos.data();
    const auto* eps = m_epsilon.data();
    const auto* tmat = m_orientation.data();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const vs::Vector cc{
            problo[0] + (i + 0.5) * dx[0],
            problo[1] + (j + 0.5) * dx[1],
            problo[2] + (k + 0.5) * dx[2],
        };

        for (int ip = 0; ip < npts; ++ip) {
            const auto dist = cc - pos[ip];
            const auto dist_local = tmat[ip] & dist;
            const auto gauss_fac = utils::gaussian3d(dist_local, eps[ip]) * vol;

            m_sum[ip] += gauss_fac;

            m_vel[ip][0] += gauss_fac * varr(i,j,k,0);
            m_vel[ip][1] += gauss_fac * varr(i,j,k,1);
            m_vel[ip][2] += gauss_fac * varr(i,j,k,2);
        }
    });

//    for (int ip = 0; ip < npts; ++ip) {
//        amrex::Print(amrex::Print::AllProcs)
//            << "[VelAvgOp] ds "
//            << dx[0] << " " << dx[1] << " " << dx[2]
//            << " pt "
//            << m_pos[ip][0] << " " << m_pos[ip][1] << " " << m_pos[ip][2]
//            << " vel "
//            << m_vel[ip][0] << " " << m_vel[ip][1] << " " << m_vel[ip][2]
//            << std::endl;
//    }
}

template <typename ActTrait>
void VelAvgOp<ActTrait>::post()
{
    auto& grid = m_data.grid();
    const int npts = m_pos.size();

    copy_to_host();

    // flatten VectorList
    amrex::Vector<amrex::Real> packed_vel(npts*3);
    for (int ip = 0, i = 0; ip < npts; ++ip) {
        for (int idim = 0; idim < 3; ++idim) {
            packed_vel[i++] = grid.vel[ip][idim];
        }
    }

//    std::set<int>::iterator procptr; 
//    amrex::Print(amrex::Print::AllProcs)
//        << "[VelAvgOp::post]"
//        << " act " << m_data.info().id
//        << " proc " << amrex::ParallelDescriptor::MyProc()
//        << " is actuator root " << m_data.info().is_root_proc
//        << std::endl;
//    if (m_data.info().is_root_proc)
//    {
//        amrex::Print(amrex::Print::AllProcs) 
//            << "turbine id=" << m_data.info().id
//            << " on procs";
//        for (procptr = m_data.info().procs.begin();
//             procptr != m_data.info().procs.end();
//             procptr++)
//            amrex::Print(amrex::Print::AllProcs) << " " << *procptr;
//        amrex::Print(amrex::Print::AllProcs) << std::endl;
//    }

    const int tag_sum = 0;
    const int tag_vel = 1;
    const size_t size = npts;
    const int nprocs = m_data.info().procs.size();
    amrex::Vector<amrex::Vector<amrex::Real>> xfer_sum(
            nprocs, amrex::Vector<amrex::Real>(npts));
    amrex::Vector<amrex::Vector<amrex::Real>> xfer_vel(
            nprocs, amrex::Vector<amrex::Real>(npts*3));
    // aggregate velocities and sum (diagnostic) onto actuator root proc
    int idx = 0;
    for (const int iproc : m_data.info().procs)
    {
        if (m_data.info().is_root_proc) // receive values
        {
            if (iproc == m_data.info().root_proc)
            {
                // don't need MPI comm on same proc
//                amrex::Print(amrex::Print::AllProcs)
//                    << "[VelAvgOp::post] act " << m_data.info().id
//                    << " proc " << amrex::ParallelDescriptor::MyProc()
//                    << " (==" << *procptr << ")"
//                    << " sets own values"
//                    << std::endl;
                xfer_sum[idx] = m_host_sum;
                xfer_vel[idx] = packed_vel;
            }
            else
            {
//                amrex::Print(amrex::Print::AllProcs)
//                    << "[VelAvgOp::post] act " << m_data.info().id
//                    << " proc " << amrex::ParallelDescriptor::MyProc()
//                    << " (==" << m_data.info().root_proc << ")"
//                    << " to receive from " << *procptr
//                    << std::endl;
                amrex::ParallelDescriptor::Recv(
                    xfer_sum[idx].data(),
                    size,
                    iproc,
                    tag_sum);
                amrex::ParallelDescriptor::Recv(
                    xfer_vel[idx].data(),
                    size*3,
                    iproc,
                    tag_vel);
            }
        }
        else if (amrex::ParallelDescriptor::MyProc() == iproc) // send
        {
//            amrex::Print(amrex::Print::AllProcs)
//                << "[VelAvgOp::post] act " << m_data.info().id
//                << " proc " << amrex::ParallelDescriptor::MyProc()
//                << " (==" << *procptr << ")"
//                << " send to " << m_data.info().root_proc
//                << std::endl;
            amrex::ParallelDescriptor::Send(
                m_host_sum.data(),
                size,
                m_data.info().root_proc,
                tag_sum);
            amrex::ParallelDescriptor::Send(
                packed_vel.data(),
                size*3,
                m_data.info().root_proc,
                tag_vel);
        }
        idx++;
    }

    if (m_data.info().is_root_proc)
    {
        // sum velocities across procs
        std::fill(packed_vel.begin(), packed_vel.end(), 0.0);
        for (int iproc = 0; iproc < nprocs; ++iproc) {
            for (int i = 0; i < npts*3; ++i) {
                packed_vel[i] += xfer_vel[iproc][i];
            }
        }

        // update diagnostic sum
        amrex::Real sum_min{9e9}, sum_max{-9e9};
        for (int ip = 0; ip < npts; ++ip) {
            m_host_sum[ip] = 0.0;
            for (int iproc = 0; iproc < nprocs; ++iproc) {
                m_host_sum[ip] += xfer_sum[iproc][ip];
            }
            if (m_host_sum[ip] < sum_min) sum_min = m_host_sum[ip];
            if (m_host_sum[ip] > sum_max) sum_max = m_host_sum[ip];
        }

        amrex::Print(amrex::Print::AllProcs)
            << "[VelAvgOp::post] actuator " << m_data.info().id
            << " integrated weight min,max "
            << sum_min << " " << sum_max << std::endl;
    }

    // scatter averaged velocities
    if (m_data.info().is_root_proc) {
        for (const int iproc : m_data.info().procs) {
            if (iproc == m_data.info().root_proc) continue;
            amrex::ParallelDescriptor::Send(
                packed_vel.data(),
                size*3,
                iproc,
                tag_vel+1);
        }
    } else {
        amrex::ParallelDescriptor::Recv(
            packed_vel.data(),
            size*3,
            m_data.info().root_proc,
            tag_vel+1);
    }

    // unpack velocities
    for (int ip = 0, i = 0; ip < npts; ++ip) {
        for (int idim = 0; idim < 3; ++idim) {
            grid.vel[ip][idim] = packed_vel[i++];
        }
    }
}

} // namespace amr_wind::actuator::ops

#endif /* VELAVGOP_H_ */
