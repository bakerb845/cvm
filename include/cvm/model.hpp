#ifndef CVM_MODEL_HPP
#define CVM_MODEL_HPP
#include <memory>
#include "cvm/enums.hpp"
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
    /// @result True indicates that the model is loaded.
    [[nodiscard]] bool isLoaded() const noexcept;

    /// @brief Destructor.
    ~Model();
    /// @brief Releases the memory and resets the class.
    void clear() noexcept;

    /// @result A pointer to the P velocity model.
    /// @throws std::rutnime_error if \c isLoaded() is false.
    [[nodiscard]] const float *getPVelocityPointer() const;
    /// @result A pointer to the S velocity model.
    /// @throws std::rutnime_error if \c isLoaded() is false.
    [[nodiscard]] const float *getSVelocityPointer() const;


    /// @result The number of grid points in x.
    /// @throws std::runtime_error if \c isLoaded() is false.
    [[nodiscard]] int getNumberOfGridPointsInX() const;
    /// @result The number of grid points in y.
    /// @throws std::runtime_error if \c isLoaded() is false.
    [[nodiscard]] int getNumberOfGridPointsInY() const;
    /// @result The number of grid points in z.
    /// @throws std::runtime_error if \c isLoaded() is false.
    [[nodiscard]] int getNumberOfGridPointsInZ() const;

    /// @result The grid spacing in x in meters.
    /// @throws std::runtime_error if \c isLoaded() is false.
    [[nodiscard]] double getGridSpacingInX() const;
    /// @result The grid spacing in y in meters.
    /// @throws std::runtime_error if \c isLoaded() is false.
    [[nodiscard]] double getGridSpacingInY() const;
    /// @result The grid spacing in z in meters.
    /// @throws std::runtime_error if \c isLoaded() is false.
    [[nodiscard]] double getGridSpacingInZ() const;


    void writeVelocities(const Options &options,
                         FileType fileType = FileType::NLL) const;
    void writeVelocities(const std::string &pFileName,
                         const std::string &sFileName,
                         FileType fileType = FileType::NLL) const;
private:
    class ModelImpl;
    std::unique_ptr<ModelImpl> pImpl;
};
}
#endif
