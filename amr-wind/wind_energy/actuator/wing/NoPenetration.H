#ifndef NO_PENETRATION_H
#define NO_PENETRATION_H

#include "amr-wind/wind_energy/actuator/wing/ActuatorWing.H"
#include "amr-wind/wind_energy/actuator/aero/AirfoilTable.H"

namespace amr_wind {
namespace actuator {

struct NoPenetrationData : public WingBaseData
{
    //int num_pts; // number of points along the span (defined in WingBaseData) 
    int num_chord_pts; // number of points along the chord (NOT IMPLEMENTED)
    RealList span_locs{0.0, 1.0};
    RealList chord_inp{1.0, 1.0};

    vs::Vector epsilon_chord;
};

struct NoPenetration : public WingType
{
    using InfoType = ActInfo;
    using GridType = ActGrid;
    using MetaType = NoPenetrationData;
    using DataType = ActDataHolder<NoPenetration>;

    static std::string identifier() { return "NoPenetration"; }
};

} // namespace actuator
} // namespace amr_wind

#endif /* NO_PENETRATION_H */
