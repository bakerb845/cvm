#ifndef CVM_SELECTION_HPP
#define CVM_SELECTION_HPP
#include <memory>
#include "cvm/enums.hpp"
namespace CVM
{
/// @brief This allows a selection (subset) of the CVM to be loaded from
///        disk instead of the entire thing. 
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Selection
{
public:
    /// @brief Default c'tor.
    Selection();
    /// @brief Copy c'tor.
    Selection(const Selection &selection);
    /// @brief Move c'tor.
    Selection(Selection &&selection) noexcept;
    /// @brief Copy assignment operator.
    /// @param[in] selection  The selection class to copy to this.
    /// @result A deep copy of the input selection.
    Selection& operator=(const Selection &selection);
    /// @brief Move assignment operator.
    /// @param[in,out] selection  The selection class whose memory will be
    ///                           moved to this.  On exit, selection's
    ///                           behavior is undefined.
    /// @result The memory from selection moved to this. 
    Selection& operator=(Selection &&selection) noexcept;
    /// @brief Destructor.
    ~Selection();
    /// @result If true then read the entire model.
    [[nodiscard]] bool readAll() const noexcept;
    /// @param[in] latlon   latlon.first is the minimum latitude in the 
    ///                     selection and latlon.second is the minimum
    ///                     longitude in the selection.
    /// @note The minimum latitude and longitude may be changed slightly
    ///       for interpolation reasons.
    /// @throws std::invalid_argument if latlon.first is less than the
    ///         minimum latitude in the CVM or latlon.second is greater
    ///         than minimum longitude in the CVM.
    void setMinimumLatitudeAndLongitude(std::pair<double, double> latlon);
    /// @param[in] latlon   latlon.first is the maximum latitude in the 
    ///                     selection and latlon.second is the maximum
    ///                     longitude in the selection.
    /// @note The maximum latitude and longitude may be changed slightly
    ///       for interpolation reasons.
    /// @throws std::invalid_argument if latlon.first is greater than the
    ///         maximum latitude in the CVM or latlon.second is greater
    ///         than the maximum longitude in the CVM.
    void setMaximumLatitudeAndLongitude(std::pair<double, double> latlon);
    /// @param[in] depths   depths.first is the minimum depth in meters
    ///                     and depths.second is the maximum depth in meters.
    /// @throw std::invalid_argument if depths.first is greater than
    ///        depths.second or depths.first or depths.second are out of
    ///        the model bounds.
    void setMinimumAndMaximumDepth(std::pair<double, double> depth);

    /// @result The starting x index when unpacking the layer.
    [[nodiscard]] int getStartPointInX(LayerIdentifier layer) const noexcept;
    /// @result The ending x index (exclusive) when unpacking the layer.
    [[nodiscard]] int getEndPointInX(LayerIdentifier layer) const noexcept;
    /// @result The starting y index when unpacking the layer.
    [[nodiscard]] int getStartPointInY(LayerIdentifier layer) const noexcept;
    /// @result The ending y index (exclusive) when unpacking the layer.
    [[nodiscard]] int getEndPointInY(LayerIdentifier layer) const noexcept;
private:
    class SelectionImpl;
    std::unique_ptr<SelectionImpl> pImpl;
};
}
#endif
