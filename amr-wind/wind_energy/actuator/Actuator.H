#ifndef ACTUATOR_H
#define ACTUATOR_H

#include "amr-wind/core/Physics.H"

/** \defgroup actuator Actuator source terms
 *  Actuator source modeling for wind turbines
 *  \ingroup wind
 */

namespace amr_wind {

class Field;

namespace actuator {

class ActuatorModel;
class ActuatorContainer;

/** Actuator line/disk modeling for wind turbines.
 *
 *  \ingroup actuator
 *
 *  This class provides an interface to model wind turbines in a wind farm as
 *  actuator line/disk source terms.
 *
 *  \sa ActuatorModel, ActuatorContainer
 */
class Actuator : public Physics::Register<Actuator>
{
public:
    static std::string identifier() { return "Actuator"; }

    explicit Actuator(CFDSim& /*sim*/);

    ~Actuator() override;

    void
    initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/) override
    {}

    void pre_init_actions() override;

    void post_init_actions() override;

    void post_regrid_actions() override;

    void pre_advance_work() override;

    void post_advance_work() override;

protected:
    //! Total number of actuator components (e.g., turbines) in the flowfield
    int num_actuators() const { return static_cast<int>(m_actuators.size()); }

    virtual void prepare_outputs();

private:
    void setup_container();

    void update_positions();

    void sample_weighted_velocities(); // if m_weighted_sampling

    void update_velocities();

    void compute_forces();

    void compute_source_term();

    void communicate_turbine_io();

    CFDSim& m_sim;

    Field& m_act_source;

    std::vector<std::unique_ptr<ActuatorModel>> m_actuators;

    std::unique_ptr<ActuatorContainer> m_container;

    //
    // ThinBody actuator options
    //
    //-control calls to ActuatorContainer object
    bool m_need_mark_faces{true};
    bool m_update_actuator_faces{false}; // at every step

    //-if true, use weighted volume sampling rather than interpolating
    // to a point with particles
    bool m_weighted_sampling{false}; // deprecated

    //-limit max source term magnitude -- for testing only, does _not_
    // account for pressure gradient contribution to RHS source
    bool m_limit_source_cfl{false};
    amrex::Real m_cfl_limit{1.0};
};

} // namespace actuator
} // namespace amr_wind

#endif /* ACTUATOR_H */
