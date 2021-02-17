#include "cvm/constants.hpp"
#include "include/pConstants.hpp"

using namespace PCVM;

Constants::Constants() = default;
Constants::~Constants() = default;

/// Grid spacing
double Constants::getGridSpacingInX(const CVM::LayerIdentifier layer) const noexcept
{
    return CVM::Constants::getGridSpacingInX(layer);
}

double Constants::getGridSpacingInY(const CVM::LayerIdentifier layer) const noexcept
{
    return CVM::Constants::getGridSpacingInY(layer);
}

double Constants::getGridSpacingInZ(const CVM::LayerIdentifier layer) const noexcept
{
    return CVM::Constants::getGridSpacingInZ(layer);
}

double Constants::getLeastCommonMultipleOfGridSpacingsInXAndY() const noexcept
{
    return CVM::Constants::getLeastCommonMultipleOfGridSpacingsInXAndY();
}

/// Grid points
int Constants::getNumberOfGridPointsInX(const CVM::LayerIdentifier layer) const noexcept
{
    return CVM::Constants::getNumberOfGridPointsInX(layer);
}

int Constants::getNumberOfGridPointsInY(const CVM::LayerIdentifier layer) const noexcept
{
    return CVM::Constants::getNumberOfGridPointsInY(layer);
}

int Constants::getNumberOfGridPointsInZ(const CVM::LayerIdentifier layer) const noexcept
{
    return CVM::Constants::getNumberOfGridPointsInZ(layer);
}

/// UTM zone
int Constants::getUTMZone() const noexcept
{
    return CVM::Constants::getUTMZone();
}

/// Number of layers
int Constants::getNumberOfLayers() const noexcept
{
    return CVM::Constants::getNumberOfLayers();
}

/// Min/max lat/lon
double Constants::getMinimumLatitude() const noexcept
{
    return CVM::Constants::getMinimumLatitude();
}

double Constants::getMinimumLongitude() const noexcept
{
    return CVM::Constants::getMinimumLongitude();
}

double Constants::getMaximumLatitude() const noexcept
{
    return CVM::Constants::getMaximumLatitude();
}

double Constants::getMaximumLongitude() const noexcept
{
    return CVM::Constants::getMaximumLongitude();
}

/// UTM origin in x and y
double Constants::getUTMOriginInX() const noexcept
{
    return CVM::Constants::getUTMOriginInX();
}

double Constants::getUTMOriginInY() const noexcept
{
    return CVM::Constants::getUTMOriginInY();
}


void PCVM::initializeConstants(pybind11::module &m)
{
    pybind11::class_<PCVM::Constants> c(m, "Constants");
    c.def(pybind11::init<> ());

    c.def("get_grid_spacing_in_x",
          &PCVM::Constants::getGridSpacingInX,
          "The grid spacing in meters in x (longitude) in a given layer in the CVM.");
    c.def("get_grid_spacing_in_y",
          &PCVM::Constants::getGridSpacingInY,
          "The grid spacing in meters in y (latitude) in a given layer in the CVM.");
    c.def("get_grid_spacing_in_z",
          &PCVM::Constants::getGridSpacingInZ,
          "The grid spacing in meters in z (depth) in a given layer in the CVM.");

    c.def("get_number_of_grid_points_in_x",
          &PCVM::Constants::getNumberOfGridPointsInX,
          "The number of grid points in x (longitude) in the CVM in a given layer.");
    c.def("get_number_of_grid_points_in_y",
          &PCVM::Constants::getNumberOfGridPointsInY,
          "The number of grid points in y (latitude) in the CVM in a given layer.");
    c.def("get_number_of_grid_points_in_z",
          &PCVM::Constants::getNumberOfGridPointsInZ,
          "The number of grid points in z (depth) in the CVM in a given layer.");

    c.def("get_utm_zone",
          &PCVM::Constants::getUTMZone,
          "The UTM zone in which the CVM is defined.");
    c.def("get_number_of_layers",
          &PCVM::Constants::getNumberOfLayers,
          "The number of layers defining the CVM.");

    c.def("get_minimum_latitude",
          &PCVM::Constants::getMinimumLatitude,
          "The smallest latitude in degrees in the CVM.");
    c.def("get_minimum_longitude",
          &PCVM::Constants::getMinimumLongitude,
          "The smallest longitude in degrees in the CVM.");

    c.def("get_maximum_latitude",
          &PCVM::Constants::getMaximumLatitude,
          "The largest latitude in degrees in the CVM.  While a layer may technically extend past the this latitude the interpolation is not well-defined.");
    c.def("get_maximum_longitude",
          &PCVM::Constants::getMaximumLongitude,
          "The largest longitude in degrees in the CVM.  While a layer may technically extend past the this longitude the interpolation is not well-defined.");

}
