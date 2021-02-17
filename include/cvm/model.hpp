#ifndef CVM_MODEL_HPP
#define CVM_MODEL_HPP
#include <memory>
namespace CVM
{
class Options;
/// @brief This class loads the CVM.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license. 
class Model
{
public:
    /// @brief Constructor.
    Model();
    /// @brief Copy constructor.
    /// @param[in] model  The model from which to initialize this class.
    Model(const Model &model);
    /// @brief Move constructor.
    /// @param[in,out] model  The model from which to initialize this class.
    ///                       On exit, model's behavior is undefined.
    Model(Model &&model) noexcept;

    /// @brief Copy assignment operator.
    /// @param[in] model  The model to copy to this.
    /// @result A deep copy of the input model. 
    Model& operator=(const Model &model);
    /// @brief Move assignment operator.
    /// @param[in,out] model  The model whose memory will be moved to this.
    ///                       On exit, model's behavior is undefined.
    /// @result The memory from model moved to this.
    Model& operator=(Model &&model) noexcept;

    /// @brief Loads the CVM.
    /// @param[in] options  The program options which contain the file names
    ///                     of all three P and S layers of the CVM.
    void load(const Options &options);

    /// @brief Destructor.
    ~Model();
    /// @brief Releases the memory and resets the class.
    void clear() noexcept;

private:
    class ModelImpl;
    std::unique_ptr<ModelImpl> pImpl;
};
}
#endif
