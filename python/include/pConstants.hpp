#ifndef PYCVM_CONSTANTS_HPP
#define PYCVM_CONSTANTS_HPP
#include <pybind11/pybind11.h>
#include "cvm/enums.hpp"
namespace PCVM
{
/// @brief Interface to the constants in the CVM library.
class Constants
{
public:
    /// C'tor
    Constants();
    /// Destructor
    ~Constants();
    /// @result The grid spacing in x in meters in the layer'th layer.
    [[nodiscard]] double getGridSpacingInX(CVM::LayerIdentifier layer) const noexcept;
    /// @result The grid spacing in y in meters in the layer'th layer.
    [[nodiscard]] double getGridSpacingInY(CVM::LayerIdentifier layer) const noexcept;
    /// @result The grid spacing in z in meters in the layer'th layer.
    [[nodiscard]] double getGridSpacingInZ(CVM::LayerIdentifier layer) const noexcept;
    /// @result Gets the least common multiple of the grid spacings
    ///         for all layers x and y in meters.
    [[nodiscard]] double getLeastCommonMultipleOfGridSpacingsInXAndY() const noexcept;
 
    /// @result The number of grid points in x in the layer'th layer. 
    [[nodiscard]] int getNumberOfGridPointsInX(CVM::LayerIdentifier layer) const noexcept;
    /// @result The number of grid points in y in the layer'th layer.
    [[nodiscard]] int getNumberOfGridPointsInY(CVM::LayerIdentifier layer) const noexcept; 
    /// @result The number of grid points in z in the layer'th layer.
    [[nodiscard]] int getNumberOfGridPointsInZ(CVM::LayerIdentifier layer) const noexcept;
   
    /// @result The UTM zone in which the CVM is defined.
    [[nodiscard]] int getUTMZone() const noexcept;

    /// @result Gets the number of layers in the CVM. 
    [[nodiscard]] int getNumberOfLayers() const noexcept;

    /// @result The smallest latitude, in degrees, in the model.
    [[nodiscard]] double getMinimumLatitude() const noexcept;
    /// @result The largest usable latitude, in degrees, in the model.
    /// @note A layer may extend beyond this latitude but interpolation
    ///       at that point is ill-defined.
    [[nodiscard]] double getMaximumLatitude() const noexcept;
    /// @result The smallest longitude, in degrees, in the model.
    [[nodiscard]] double getMinimumLongitude() const noexcept;
    /// @result The largest usable longitude, in degrees, in the model.
    /// @note A layer may extend beyond this longitude but interpolation
    ///       at that point is ill-defined.
    [[nodiscard]] double getMaximumLongitude() const noexcept;


    /// @result The x UTM (zone 10) of the lower left corner of the model
    ///         in meters.
    [[nodiscard]] double getUTMOriginInX() const noexcept; 
    /// @result The y UTM (zone 10) of the lower left corner of the model
    ///         in meters.
    [[nodiscard]] double getUTMOriginInY() const noexcept;
};
void initializeConstants(pybind11::module &m);
}
#endif

