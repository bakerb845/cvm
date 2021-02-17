#ifndef CVM_CONSTANTS_HPP
#define CVM_CONSTANTS_HPP
#include <memory>
#include "cvm/enums.hpp"
namespace CVM
{
/// @brief Defines constant grid parameters defining the CVM.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Constants
{
public:
    /// @result The grid spacing in x in meters in the layer'th layer.
    [[nodiscard]] static double getGridSpacingInX(LayerIdentifier layer) noexcept;
    /// @result The grid spacing in y in meters in the layer'th layer.
    [[nodiscard]] static double getGridSpacingInY(LayerIdentifier layer) noexcept;
    /// @result The grid spacing in z in meters in the layer'th layer.
    [[nodiscard]] static double getGridSpacingInZ(LayerIdentifier layer) noexcept;
    /// @result Gets the least common multiple of the grid spacings
    ///         for all layers x and y in meters.
    [[nodiscard]] static double getLeastCommonMultipleOfGridSpacingsInXAndY() noexcept;
 
    /// @result The number of grid points in x in the layer'th layer. 
    [[nodiscard]] static int getNumberOfGridPointsInX(LayerIdentifier layer) noexcept;
    /// @result The number of grid points in y in the layer'th layer.
    [[nodiscard]] static int getNumberOfGridPointsInY(LayerIdentifier layer) noexcept; 
    /// @result The number of grid points in z in the layer'th layer.
    [[nodiscard]] static int getNumberOfGridPointsInZ(LayerIdentifier layer) noexcept;

    /// @result The maximum usable number of grid points in x in the layer'th layer.
    /// @bug This is slightly inconsistent with \c getMaximumLongitude() since 
    ///      \c getMaximumLongitude() produces a slightly larger max longitude.
    [[nodiscard]] static int getMaximumUsableNumberOfGridPointsInX(LayerIdentifier layer) noexcept;
    /// @result The maximum usable number of grid points in y in the layer'th layer.
    /// @bug This is slightly inconsistent with \c getMaximumLatitude() since
    ///      \c getMaximumLatitude() produces a slightly larger max latitude.
    [[nodiscard]] static int getMaximumUsableNumberOfGridPointsInY(LayerIdentifier layer) noexcept;
   
    /// @result The UTM zone in which the CVM is defined.
    [[nodiscard]] static int getUTMZone() noexcept;

    /// @result Gets the number of layers in the CVM. 
    [[nodiscard]] static int getNumberOfLayers() noexcept;

    /// @result The smallest latitude, in degrees, in the model.
    [[nodiscard]] static double getMinimumLatitude() noexcept;
    /// @result The largest usable latitude, in degrees, in the model.
    /// @note A layer may extend beyond this latitude but interpolation
    ///       at that point is ill-defined.
    [[nodiscard]] static double getMaximumLatitude() noexcept;
    /// @result The smallest longitude, in degrees, in the model.
    [[nodiscard]] static double getMinimumLongitude() noexcept;
    /// @result The largest usable longitude, in degrees, in the model.
    /// @note A layer may extend beyond this longitude but interpolation
    ///       at that point is ill-defined.
    [[nodiscard]] static double getMaximumLongitude() noexcept;


    /// @result The x UTM (zone 10) of the lower left corner of the model
    ///         in meters.
    [[nodiscard]] static double getUTMOriginInX() noexcept; 
    /// @result The y UTM (zone 10) of the lower left corner of the model
    ///         in meters.
    [[nodiscard]] static double getUTMOriginInY() noexcept;

    /// @result The top-most depth of the layer in meters.
    [[nodiscard]] static double getLayerStartDepth(LayerIdentifier layer) noexcept;
    /// @result The bottom-most depth of the layer in meters. 
    [[nodiscard]] static double getLayerEndDepth(LayerIdentifier layer) noexcept;
};
}
#endif
