#ifndef CVM_GEODETIC_HPP
#define CVM_GEODETIC_HPP
#include <utility>
namespace CVM::Geodetic
{
/// @param[in] latlon  latlon.first is the latitude in degrees and 
///                    latlon.second is the longitude in degrees.
/// @param[in] utmProjectionZone  The UTM projection zone.  The default will
///                               be appropriate for the CVM.
/// @result result.first is the UTM in x in meters and result.second is the
///         UTM in y in meters corresponding to the given lat/lon. 
std::pair<double, double> latitudeLongitudeToUTM(
     const std::pair<double, double> &latlon,
     int utmProjectionZone = 10);
/// @param[in] xy   xy.first is the UTM x position in meters 
///                 and xy.second is the UTM y position in meters.
/// @param[in] utmProjectionZone  The UTM projection zone.  The default will
///                               be appropriate for the CVM.
/// @result result.first is the latitude in degrees and result.second
///         is the longitude in degrees corresponding to the given UTM.
std::pair<double, double> utmToLatitudeLongitude(
     const std::pair<double, double> &xy,
     int utmProjectionZone = 10);
}
#endif
