#ifndef PYCVM_GEODETIC_HPP
#define PYCVM_GEODETIC_HPP
#include <utility>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace PCVM
{
/// @brief Interface to the simple geodetic functionality in the CVM library.
class Geodetic
{
public:
    /// C'tor for zone 10
    Geodetic(int zone = 10);
    Geodetic(const Geodetic &geodetic) = default;
    Geodetic(Geodetic &&geodetic) noexcept = default;
    Geodetic& operator=(const Geodetic &geodetic) = default;
    Geodetic& operator=(Geodetic &&geodetic) noexcept = default;
    /// Destructor
    ~Geodetic();
    /// Lat/lon to UTM 
    std::pair<double, double> latitudeLongitudeToUTM(const std::pair<double, double> &latlon);
    /// UTM to lat/lon
    std::pair<double, double> utmToLatitudeLongitude(const std::pair<double, double> &utm);
private:
    int mZone = 10;
};
void initializeGeodetic(pybind11::module &m);
}
#endif
