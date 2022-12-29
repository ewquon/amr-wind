#include "amr-wind/wind_energy/actuator/wing/ThinBody.H"
#include "amr-wind/wind_energy/actuator/wing/thin_body_ops.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"

namespace amr_wind::actuator {

template class ActModel<ThinBody, ActSrcLine>;

} // namespace amr_wind::actuator
