#ifndef THIN_BODY_OPS_H
#define THIN_BODY_OPS_H

#include "amr-wind/wind_energy/actuator/wing/ThinBody.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "amr-wind/wind_energy/actuator/FLLCOp.H"

namespace amr_wind::actuator {
namespace wing {

/** Parse common data for modeling fixed wings as actuators

    Copied from wing_ops.H
    ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
 */
void read_inputs(
    WingBaseData& /*wdata*/, ActInfo& /*info*/, const utils::ActParser& /*pp*/);

/** Initialize core data structures when modeling fixed wings
 */
void init_data_structures(WingBaseData& /*wdata*/, ActGrid& /*grid*/);
void init_data_structures_no_edge(WingBaseData& /*wdata*/, ActGrid& /*grid*/);

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

template <typename T>
ComponentView make_component_view(typename T::DataType& data)
{
    auto& grid = data.grid();
    auto& meta = data.meta();
    ComponentView view;
    view.pos = ::amr_wind::utils::slice(grid.pos, 0, meta.num_pts);
    view.vel_pos = ::amr_wind::utils::slice(grid.vel_pos, 0, meta.num_pts);
    view.force = ::amr_wind::utils::slice(grid.force, 0, meta.num_pts);
    view.epsilon = ::amr_wind::utils::slice(grid.epsilon, 0, meta.num_pts);
    view.orientation =
        ::amr_wind::utils::slice(grid.orientation, 0, meta.num_pts);
    view.chord = ::amr_wind::utils::slice(meta.chord, 0, meta.num_pts);
    view.vel_rel = ::amr_wind::utils::slice(meta.vel_rel, 0, meta.num_pts);
    view.vel = ::amr_wind::utils::slice(grid.vel, 0, meta.num_pts);

    return view;
}

/** based on wing_ops.cpp
    ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓

  - actuator points do not extend to edge; rationale is that the body force
    distribution should be able to represent a pressure jump over the thin 
    body that vanishes at the edges, so the actuator forces at the end points
    should tend to 0
  - all actuator points have the same dx
 */
void init_data_structures_no_edge(WingBaseData& wdata, ActGrid& grid)
{
    int npts = wdata.num_pts;
    grid.resize(npts);

    // Wing span
    auto wspan = wdata.end - wdata.start;
    // Compute transformation matrix as a quaternion rotation about wing span
    auto tmat = vs::quaternion(wspan, wdata.pitch);
    // Epsilon (chord, thickness, span) -> (chord, span, thickness)
    const auto& epsin = wdata.eps_inp;
    vs::Vector eps{epsin.x(), epsin.z(), epsin.y()};

    // Equal spacing along span
    auto dx = (1.0 / static_cast<amrex::Real>(npts)) * wspan;

    for (int i = 0; i < npts; ++i) {
        grid.pos[i] = wdata.start + (static_cast<amrex::Real>(i)+0.5) * dx;
    }

    // Initialize remaining data
    grid.epsilon.assign(npts, eps);
    grid.orientation.assign(npts, tmat);
    grid.force.assign(npts, vs::Vector::zero());
    grid.vel_pos.assign(grid.pos.begin(), grid.pos.end());
    grid.vel.assign(npts, vs::Vector::zero());

    wdata.dx.assign(npts, vs::mag(dx));
    wdata.vel_rel.assign(npts, vs::Vector::zero());
    wdata.aoa.assign(npts, 0.0);
    wdata.cl.assign(npts, 0.0);
    wdata.cd.assign(npts, 0.0);
}

} // namespace wing

namespace ops {

/** based on fixed_wing_ops.H
    ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓

  - removed airfoil table lookup
  - add params `gain`, `num_chord_points` (not yet implemented)
  - calculate max allowable point force to impose CFL limit
 */

template <>
struct ReadInputsOp<ThinBody, ActSrcLine>
{
    void operator()(ThinBody::DataType& data, const utils::ActParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();
        pp.get("num_span_points", wdata.num_pts);
        // pp.get("num_chord_points", wdata.num_chord_pts);
        pp.query("actuators_at_nodes", wdata.act_nodes);
        pp.get("start", wdata.start);
        pp.get("end", wdata.end);
        pp.query("gain", wdata.gain);
        pp.query("epsilon", wdata.eps_inp);
        pp.query("epsilon_chord", wdata.epsilon_chord);
        pp.query("pitch", wdata.pitch);
        pp.queryarr("span_locs", wdata.span_locs);
        pp.queryarr("chord", wdata.chord_inp);

        bool use_fllc = false;
        pp.query("fllc", use_fllc);
        if (use_fllc) {
            wdata.fllc = std::make_unique<FLLCData>();
            FLLCParse(pp, *(wdata.fllc));
        }

        if (!pp.contains("epsilon") && !pp.contains("epsilon_chord")) {
            amrex::Abort(
                "Actuator thin body requires specification of one or both "
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
            amrex::max(max_eps, max_epsc) * max_chord * 50.0; //3.0; -- ad hoc value to produe identical result as serial run
        //amrex::Print() << "search_radius = " << search_radius << std::endl;
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
struct InitDataOp<ThinBody, ActSrcLine>
{
    void operator()(ThinBody::DataType& data)
    {
        // this sets up the actuator grid
        if (data.meta().act_nodes) {
            // force points extend to edge of actuator grid; first and
            // last points have half the segment width (actuator line
            // model default)
            wing::init_data_structures(data.meta(), data.grid());
        } else {
            // force pts are at the center of each actuator segment
            wing::init_data_structures_no_edge(data.meta(), data.grid());
        }

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

        // Create field to use in CV analysis adjacent to the actuator surface
        // ncomp = 3, ngrow = 1
        data.sim().repo().declare_field("mom_sum",3,1);

        // Create masks used to zero the surface velocity/flux
        // ncomp = 1, ngrow = 0
        data.sim().repo().declare_int_field("mom_xface_mask",1,0,1,amr_wind::FieldLoc::XFACE);
        data.sim().repo().declare_int_field("mom_yface_mask",1,0,1,amr_wind::FieldLoc::YFACE);
        data.sim().repo().declare_int_field("mom_zface_mask",1,0,1,amr_wind::FieldLoc::ZFACE);

        // Create field to hold the normal velocity corrections for cells
        // adjacent to the surface.
        data.sim().repo().declare_field("norm_vel_corr",3,0,1,amr_wind::FieldLoc::CELL);

        // Create masks for overset-like solves
        // TODO: Make sure there's no conflict with existing fields in overset
        // infrastructure
        int nghost = data.sim().pde_manager().num_ghost_state();
        data.sim().repo().declare_int_field("mask_cell",1,nghost);
        data.sim().repo().declare_int_field("mask_node",1,nghost,1,FieldLoc::NODE);
    }
};

/** Copied from wing_ops.H
    ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
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
    void operator()(typename ActTrait::DataType& data)
    {
        auto& meta = data.meta();
        if (meta.fllc) {
            FLLCOp()(meta.component_view, *(meta.fllc));
        }
    }
};


/** specific to ThinBody actuator
    ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓

    Actuator::post_init_actions, Actuator::pre_advance_work
    --> Actuator::compute_forces
        --> ActModel::compute_forces -- if data.info().actuator_in_proc
            --> ops::ComputeForceOp()
            --> ops::ActSrcOp.setup_op() -- copies newly computed forces to device
 */

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

        //amrex::Real dt = data.sim().time().deltaT();

        amrex::Real total_lift = 0.0;
        amrex::Real total_drag = 0.0;
        for (int ip = 0; ip < npts; ++ip) {
            const auto& tmat = grid.orientation[ip];
            // Effective velocity at the wing control point in local frame
            auto wvel = tmat & grid.vel[ip];
            // Set spanwise component to zero to get a pure 2D velocity
            wvel.y() = 0.0;

            if (DEBUG && data.info().is_root_proc) {
                //if ((data.info().id == 0) && (ip == 0))
                amrex::Print(amrex::Print::AllProcs)
                    << "[ComputeForceOp] actuator " << data.info().id
                    << " pt "<<ip<<" at ("<<grid.pos[ip].x()<<","<<grid.pos[ip].y()<<","<<grid.pos[ip].z()<< ")"
                    << " wvel=("<<wvel.x()<<","<<wvel.y()<<","<<wvel.z()<<")"
                    << " force=("<<grid.force[ip].x()<<","<<grid.force[ip].y()<<","<<grid.force[ip].z()<<")"
                    << std::endl;
            }

            // Compute force on fluid to achieve no penetration
            // assuming unit density: f ~ dv/dt
            auto wvel_normal = wvel;
            wvel_normal.x() = 0.0;

            // grid.force is set to the CV analysis estimate with ActModel::update_fields()
            //grid.force[ip] = wdata.gain * -(wvel_normal & tmat) / dt; // in inertial frame

            // Compute 2-D lift, drag assuming unit density
            amrex::Real cl, cd;
            const auto vmag = vs::mag(wvel);
            const auto aoa = std::atan2(wvel.z(), wvel.x());
            const auto qval = 0.5 * vmag * vmag * chord[ip] * dx[ip];
            const auto drag_dir = wvel.unit();
            const auto lift_dir = drag_dir ^ tmat.y();
            amrex::Real lift = grid.force[ip] & lift_dir;
            amrex::Real drag = grid.force[ip] & drag_dir;
            cl = lift / qval;
            cd = drag / qval;

            // Assign values for output
            wdata.vel_rel[ip] = wvel;
            wdata.aoa[ip] = amr_wind::utils::degrees(aoa);
            wdata.cl[ip] = cl;
            wdata.cd[ip] = cd;

            total_lift += lift;
            total_drag += drag;

            // Body force should _not_ be needed with new interface treatment
            // Note: body force in actuator output will need to be
            // reconstructed from cl, cd, vrel, and aoa
            grid.force[ip] = amr_wind::vs::Vector::zero();
        }

        wdata.lift = total_lift;
        wdata.drag = total_drag;

        if (wdata.fllc) {
            if (!(wdata.fllc->initialized)) {
                wdata.component_view =
                    amr_wind::actuator::wing::make_component_view<ActTrait>(
                        data);
                FLLCInit(
                    *(wdata.fllc.get()), wdata.component_view,
                    wdata.epsilon_chord[0]);
            }
        }
    }
};

template <typename ActTrait>
struct ScaleForceOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if<
        std::is_base_of<WingType, ActTrait>::value &&
        UseDefaultOp<ActTrait, ActSrcLine>::compute_force>::type>
{
    void operator()(typename ActTrait::DataType& data, amrex::Real gain)
    {
        auto& grid = data.grid();
        auto& wdata = data.meta();
        const int npts = wdata.num_pts;
        for (int ip = 0; ip < npts; ++ip) {
            grid.force[ip] *= gain;
            wdata.cl[ip] *= gain;
            wdata.cd[ip] *= gain;
        }
        wdata.lift *= gain;
        wdata.drag *= gain;
    }
};

/** Copied from wing_ops.H
    ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
 */

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
} // namespace amr_wind::actuator

#endif /* THIN_BODY_OPS_H */
