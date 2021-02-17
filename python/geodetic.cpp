#include "cvm/geodetic.hpp"
#include "include/pGeodetic.hpp"

using namespace PCVM;

Geodetic::Geodetic(const int zone) :
    mZone(zone)
{
}

Geodetic::~Geodetic() = default;

std::pair<double, double> Geodetic::latitudeLongitudeToUTM(
    const std::pair<double, double> &latLon)
{
    return CVM::Geodetic::latitudeLongitudeToUTM(latLon, mZone);
}

std::pair<double, double> Geodetic::utmToLatitudeLongitude(
    const std::pair<double, double> &utm)
{
     return CVM::Geodetic::utmToLatitudeLongitude(utm, mZone);
}

void PCVM::initializeGeodetic(pybind11::module &m)
{
    pybind11::class_<PCVM::Geodetic> g(m, "Geodetic");
    g.def(pybind11::init<> ());
    g.def("latitude_longitude_to_utm",
          &PCVM::Geodetic::latitudeLongitudeToUTM,
          "Converts latitude and longitude in degrees to a UTM in meters in UTM Zone 10.");
    g.def("utm_to_latitude_longitude",
          &PCVM::Geodetic::utmToLatitudeLongitude,
          "Converts a UTM to a latitude and longitude in degrees");
}
