#ifndef NO_PENETRATION_OPS_H
#define NO_PENETRATION_OPS_H

#include "amr-wind/wind_energy/actuator/wing/NoPenetration.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

namespace amr_wind {
namespace actuator {
namespace wing {

/** Parse common data for modeling fixed wings as actuators

    Copied from wing_ops.H
    vvvvvvvvvvvvvvvvvvvvvv
 */
void read_inputs(
    WingBaseData& /*wdata*/, ActInfo& /*info*/, const utils::ActParser& /*pp*/);

/** Initialize core data structures when modeling fixed wings
 */
void init_data_structures(WingBaseData& /*wdata*/, ActGrid& /*grid*/);

void prepare_netcdf_file(
    const std::string& /*ncfile*/,
    const WingBaseData& /*meta*/,
    const ActInfo& /*info*/,
    const ActGrid& /*grid*/);

void write_netcdf(
    const std::string& /*ncfile*/,
    const WingBaseData& /*meta*/,
    const ActInfo& /*info*/,
    const ActGrid& /*grid*/,
    const amrex::Real /*time*/);

} // namespace wing


/** NoPenetration actuator specific
    vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 */
namespace ops {

template <>
struct ReadInputsOp<NoPenetration, ActSrcLine>
{
    void operator()(NoPenetration::DataType& data, const utils::ActParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();
        pp.get("num_span_points", wdata.num_pts);
        pp.get("num_chord_points", wdata.num_chord_pts);
        pp.get("start", wdata.start);
        pp.get("end", wdata.end);
        pp.query("gain", wdata.gain);
        pp.query("epsilon", wdata.eps_inp);
        pp.query("epsilon_chord", wdata.epsilon_chord);
        pp.get("pitch", wdata.pitch);
        pp.queryarr("span_locs", wdata.span_locs);
        pp.queryarr("chord", wdata.chord_inp);

        if (!pp.contains("epsilon") && !pp.contains("epsilon_chord")) {
            amrex::Abort(
                "Actuator fixed wing requires specification of one or both "
                "of 'epsilon' or 'epsilon_chord'");
        }

        AMREX_ALWAYS_ASSERT(wdata.span_locs.size() == wdata.chord_inp.size());

        amrex::Real max_chord =
            *std::max_element(wdata.chord_inp.begin(), wdata.chord_inp.end());
        amrex::Real max_eps =
            *std::max_element(wdata.eps_inp.begin(), wdata.eps_inp.end());
        amrex::Real max_epsc = *std::max_element(
            wdata.epsilon_chord.begin(), wdata.epsilon_chord.end());
        amrex::Real search_radius =
            amrex::max(max_eps, max_epsc) * max_chord * 3.0;
        const auto& p1 = wdata.start;
        const auto& p2 = wdata.end;
        // clang-format off
        info.bound_box = amrex::RealBox(
            amrex::min(p1.x(), p2.x()) - search_radius,
            amrex::min(p1.y(), p2.y()) - search_radius,
            amrex::min(p1.z(), p2.z()) - search_radius,
            amrex::max(p1.x(), p2.x()) + search_radius,
            amrex::max(p1.y(), p2.y()) + search_radius,
            amrex::max(p1.z(), p2.z()) + search_radius
        );
        // clang-format on
    }
};

template <>
struct InitDataOp<NoPenetration, ActSrcLine>
{
    void operator()(NoPenetration::DataType& data)
    {
        wing::init_data_structures(data.meta(), data.grid());

        auto& meta = data.meta();
        {
            const int npts = data.meta().num_pts;
            auto& grid = data.grid();
            // Process chord information
            const auto wlen = vs::mag(grid.pos.back() - grid.pos.front());
            RealList wx(npts);
            for (int i = 0; i < npts; ++i) {
                wx[i] = vs::mag(grid.pos[i]) / wlen;
            }
            meta.chord.resize(npts);
            ::amr_wind::interp::linear_monotonic(
                meta.span_locs, meta.chord_inp, wx, meta.chord);
            // clang-format off
            meta.epsilon_chord = {
                meta.epsilon_chord.x(), meta.epsilon_chord.z(),
                meta.epsilon_chord.y()};
            // clang-format on
            for (int i = 0; i < npts; ++i) {
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    const auto eps = meta.epsilon_chord[n] * meta.chord[i];
                    grid.epsilon[i][n] = amrex::max(meta.eps_inp[n], eps);
                }
            }
        }
    }
};

/* Copied from wing_ops.H -- removed airfoil lookup from ComputeForceOp
   vvvvvvvvvvvvvvvvvvvvvv
*/

template <typename ActTrait, typename SrcTrait>
struct UseDefaultOp
{
    static constexpr bool update_pos = true;
    static constexpr bool update_vel = true;
    static constexpr bool compute_force = true;
    static constexpr bool process_outputs = true;
};

template <typename ActTrait>
struct UpdatePosOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if<
        std::is_base_of<WingType, ActTrait>::value &&
        UseDefaultOp<ActTrait, ActSrcLine>::update_pos>::type>
{
    void operator()(typename ActTrait::DataType& /*unused*/) {}
};

template <typename ActTrait>
struct UpdateVelOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if<
        std::is_base_of<WingType, ActTrait>::value &&
        UseDefaultOp<ActTrait, ActSrcLine>::update_vel>::type>
{
    void operator()(typename ActTrait::DataType& /*unused*/) {}
};

template <typename ActTrait>
struct ComputeForceOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if<
        std::is_base_of<WingType, ActTrait>::value &&
        UseDefaultOp<ActTrait, ActSrcLine>::compute_force>::type>
{
    void operator()(typename ActTrait::DataType& data)
    {
        auto& grid = data.grid();
        auto& wdata = data.meta();
        const int npts = wdata.num_pts;
        const auto& dx = wdata.dx;
        const auto& chord = wdata.chord;

        const bool DEBUG = false;

        amrex::Real dt = data.sim().time().deltaT();

        amrex::Real total_lift = 0.0;
        amrex::Real total_drag = 0.0;
        for (int ip = 0; ip < npts; ++ip) {
            const auto& tmat = grid.orientation[ip];
            // Effective velocity at the wing control point in local frame
            auto wvel = tmat & grid.vel[ip];
            // Set spanwise component to zero to get a pure 2D velocity
            wvel.y() = 0.0;

            const auto vmag = vs::mag(wvel);
            const auto aoa = std::atan2(wvel.z(), wvel.x()); // not true angle of attack
            amrex::Real cl, cd;

            // Compute force on fluid to achieve no penetration
            // assuming unit density
            auto wvel_normal = wvel;
            wvel_normal.x() = 0.0;
            wvel_normal.y() = 0.0;
            grid.force[ip] = wdata.gain * -(wvel_normal & tmat) / dt;

            if (DEBUG) {
                amrex::Print(amrex::Print::AllProcs) << "act pt " << ip << ":"
                    << " grid.ori=" << grid.orientation[ip]
                    << " grid.vel=" << grid.vel[ip]
                    //<< " lift_dir=" << lift_dir
                    //<< " drag_dir=" << drag_dir
                    << " wvel=" << wvel
                    << " force=" << grid.force[ip]
                    << std::endl;
            }

            // Determine unit vector parallel and perpendicular to velocity vector
            const auto drag_dir = wvel.unit() & tmat;
            const auto lift_dir = drag_dir ^ tmat.y();

            // Compute lift, drag assuming unit density
            const auto qval = 0.5 * vmag * vmag * chord[ip] * dx[ip];
            const auto lift = 0.0; //qval * cl;
            const auto drag = 0.0; //qval * cd;
            cl = lift / qval;
            cd = drag / qval;

            // Assign values for output
            wdata.vel_rel[ip] = wvel;
            wdata.aoa[ip] = amr_wind::utils::degrees(aoa);
            wdata.cl[ip] = cl;
            wdata.cd[ip] = cd;

            total_lift += lift;
            total_drag += drag;
        }

        wdata.lift = total_lift;
        wdata.drag = total_drag;
    }
};

template <typename ActTrait, typename SrcTrait>
struct ProcessOutputsOp<
    ActTrait,
    SrcTrait,
    typename std::enable_if<
        std::is_base_of<WingType, ActTrait>::value &&
        UseDefaultOp<ActTrait, ActSrcLine>::process_outputs>::type>
{
private:
    typename ActTrait::DataType& m_data;

    std::string m_out_dir;
    std::string m_nc_filename;
    int m_out_freq{10};

public:
    explicit ProcessOutputsOp(typename ActTrait::DataType& data) : m_data(data)
    {}

    void read_io_options(const utils::ActParser& pp)
    {
        pp.query("output_frequency", m_out_freq);
    }

    void prepare_outputs(const std::string& out_dir)
    {
        m_nc_filename = out_dir + "/" + m_data.info().label + ".nc";
        wing::prepare_netcdf_file(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid());
    }

    void write_outputs()
    {
        const auto& time = m_data.sim().time();
        const int tidx = time.time_index();
        if (tidx % m_out_freq != 0) {
            return;
        }

        wing::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid(),
            time.new_time());
    }
};

} // namespace ops
} // namespace actuator
} // namespace amr_wind

#endif /* NO_PENETRATION_OPS_H */
