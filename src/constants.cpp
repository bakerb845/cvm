#include <iostream>
#include <array>
#include "cvm/constants.hpp"

using namespace CVM;

namespace
{

int layerToInteger(const LayerIdentifier layer)
{
    return static_cast<int> (layer);
}

}

// Gets the number of layers in the CVM
int Constants::getNumberOfLayers() noexcept
{
    return 3;
}

// Grid spacing in x, y, and z in each layer
double Constants::getGridSpacingInX(const LayerIdentifier layer) noexcept
{
    std::array<double, 3> dx{200, 300, 900};
    auto idx = layerToInteger(layer);
    return dx[idx];
}

double Constants::getGridSpacingInY(const LayerIdentifier layer) noexcept
{
    std::array<double, 3> dy{200, 300, 900};
    auto idx = layerToInteger(layer);
    return dy[idx];
}

double Constants::getGridSpacingInZ(const LayerIdentifier layer) noexcept
{
    std::array<double, 3> dz{100, 300, 900};
    auto idx = layerToInteger(layer);
    return dz[idx];
}

// Grid points in x, y, and z in each layer
int Constants::getNumberOfGridPointsInX(const LayerIdentifier layer) noexcept
{
    std::array<int, 3> nx{3271, 2181,  727};
    auto idx = layerToInteger(layer);
    return nx[idx];
} 

int Constants::getNumberOfGridPointsInY(const LayerIdentifier layer) noexcept
{
    std::array<int, 3> ny{5367, 3578, 1193};
    auto idx = layerToInteger(layer);
    return ny[idx];
} 

int Constants::getMaximumUsableNumberOfGridPointsInX(
    const LayerIdentifier layer) noexcept
{
    std::array<int, 3> nx{3267, 2178, 726};
    auto idx = layerToInteger(layer);
    return nx[idx];
}

int Constants::getMaximumUsableNumberOfGridPointsInY(
    const LayerIdentifier layer) noexcept
{
    std::array<int, 3> ny{5247, 3498, 1166};
    auto idx = layerToInteger(layer);
    return ny[idx];
}

int Constants::getNumberOfGridPointsInZ(const LayerIdentifier layer) noexcept
{
    std::array<int, 3> nz{  13,   29,   55};
    auto idx = layerToInteger(layer);
    return nz[idx];
}

// UTM zone
int Constants::getUTMZone() noexcept
{
    return 10;
}

double Constants::getMinimumLatitude() noexcept
{
    return 40.2;
}

double Constants::getMinimumLongitude() noexcept
{
    return -129; //231.0;
}

double Constants::getMinimumDepth() noexcept
{
    return 0;
}

double Constants::getMaximumDepth() noexcept
{
    return 59400; 
}

double Constants::getLeastCommonMultipleOfGridSpacingsInXAndY() noexcept
{
    return 1800;
}

double Constants::getMaximumLongitude() noexcept
{
    return -121.01033362912;
}

double Constants::getMaximumLatitude() noexcept
{
    return 49.796693895676;
}

double Constants::getUTMOriginInX() noexcept
{
    return -10800.0;
}

double Constants::getUTMOriginInY() noexcept
{
    return 4467300.0;
}

double Constants::getLayerStartDepth(const LayerIdentifier layer) noexcept
{
    auto il = layerToInteger(layer);
    std::array<double, 3> depths{0, 1500, 10800};
    return depths[il];
}

double Constants::getLayerEndDepth(const LayerIdentifier layer) noexcept
{
    auto il = layerToInteger(layer);
    std::array<double, 3> depths{1200, 9900, 59400};
    return depths[il];
}
