#ifndef CVM_ENUMS_HPP
#define CVM_ENUMS_HPP
namespace CVM
{
/// @brief The CVM was traditionally stored in three separated files.
///        A fine, top layer for near-surface heterogeneities.  
///        A coarser layer for intermediate crustal parameters.
///        And a very coarse lower-crust/mantle layer.
enum class LayerIdentifier
{
    TOP = 0,    /// The upper-crust layer.
    MIDDLE = 1, /// The middle-crust layer.
    BOTTOM = 2  /// The lower-crust/mantle layer.
};
}
#endif
