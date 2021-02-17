#ifndef CVM_LAYER_HPP
#define CVM_LAYER_HPP
#include <string>
#include <memory>
#include "cvm/enums.hpp"
namespace CVM
{
class Selection;
/// @brief Reads a CVM model.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
template<LayerIdentifier E>
class Layer
{
public:
    /// @brief Constructor.
    Layer();
    /// @brief Destructor.
    ~Layer();
    /// @brief Releases memory and resets the layer.
    void clear() noexcept;
    /// @brief Loads the CVM layer from disk.
    /// @param[in] fileName   The name of the file containing the layer.
    /// @param[in] selection  Defines the subset of the velocity model to read.
    /// @param[in] isP        True indicates this is a P-velocity layer.
    ///                       False indicates this is an S-velocity layer. 
    /// @throws std::invalid_argument if the file does not exist.
    /// @throws std::runtime_error if the file is not appropriately formatted.
    void load(const std::string &fileName, const Selection  &selection, const bool isP);
    /// @result A pointer to the P velocity (km/s) for this layer.
    /// @throw std::invalid_argument if \c havePVelocity() is false. 
    [[nodiscard]] const float *getPVelocityPointer() const;
    /// @result A pointer to the S velocity (km/s) for this layer.
    /// @throw std::invalid_argument if \c haveSVelocity() is false. 
    [[nodiscard]] const float *getSVelocityPointer() const;
    /// @result A pointer to the density (gm/cm**3) for this layer.
    /// @throw std::invalid_argument if \c havePVelocity() is false. 
    [[nodiscard]] const float *getDensityPointer() const;
    /// @result True indicates the P velocity is loaded.
    [[nodiscard]] bool havePVelocity() const noexcept;
    /// @result True indicates the S velocity is loaded.
    [[nodiscard]] bool haveSVelocity() const noexcept;
    /// @result The layer identifier.
    [[nodiscard]] LayerIdentifier getLayerIdentifier() const noexcept;
    /// @result The number of grid points in x.
    /// @throws std::invalid_argument if a model has not been loaded. 
    [[nodiscard]] int getNumberOfGridPointsInX() const;
    /// @result The number of grid points in y.
    /// @throws std::invalid_argument if a model has not been loaded.
    [[nodiscard]] int getNumberOfGridPointsInY() const;
    /// @result The number of grid points in z.
    /// @throws std::invalid_argument if a model has not been loaded.
    [[nodiscard]] int getNumberOfGridPointsInZ() const;

    /// @result result.first is the UTM x position in meters of the lower
    ///         left corner of the model and result.second is the UTM
    ///         y position in meters of the lower left corner of the
    ///         model.
    /// @throws std::invalid_argument if a model has not been loaded.
    [[nodiscard]] std::pair<double, double> getSouthWestCornerUTM() const;
    /// @result result.first is the UTM x position in meters of the upper 
    ///         right corner of the model and result.second is the UTM
    ///         y position in meters of the upper right corner of the
    ///         model.
    /// @throws std::invalid_argument if a model has not been loaded.
    [[nodiscard]] std::pair<double, double> getNorthEastCornerUTM() const;
private:
    class LayerImpl;
    std::unique_ptr<LayerImpl> pImpl;
};
}
#endif
