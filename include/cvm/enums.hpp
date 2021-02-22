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
/// @brief The CVM can be written in a few foramts that are convenient for
///        different things.  For example, legacy VTK and ParaView can be
///        used to visualize the model whereas NLL can be used for subseqeuent
///        locations.
enum class FileType
{
    NLL = 0,    /// NonLinLoc file.
    VTK = 1     /// VTK legacy file format.
};
}
#endif
