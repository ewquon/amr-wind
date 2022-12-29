#ifndef NO_THINBODY_H
#define NO_THINBODY_H

#include "amr-wind/wind_energy/actuator/wing/ActuatorWing.H"

namespace amr_wind {
namespace actuator {

struct ThinBodyData : public WingBaseData
{
    //int num_pts; // number of points along the span (defined in WingBaseData) 
    int num_chord_pts; // number of points along the chord (NOT IMPLEMENTED)
    RealList span_locs{0.0, 1.0};
    RealList chord_inp{1.0, 1.0};

    amrex::Real gain = 1.0;

    vs::Vector epsilon_chord;
};

struct ThinBody : public WingType
{
    using InfoType = ActInfo;
    using GridType = ActGrid;
    using MetaType = ThinBodyData;
    using DataType = ActDataHolder<ThinBody>;

    static std::string identifier() { return "ThinBody"; }
};

} // namespace actuator
} // namespace amr_wind

#endif /* NO_THINBODY_H */
