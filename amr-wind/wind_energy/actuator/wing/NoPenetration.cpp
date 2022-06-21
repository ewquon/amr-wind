#include "amr-wind/wind_energy/actuator/wing/NoPenetration.H"
#include "amr-wind/wind_energy/actuator/wing/no_penetration_ops.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"

namespace amr_wind {
namespace actuator {

template class ActModel<NoPenetration, ActSrcLine>;

} // namespace actuator
} // namespace amr_wind
