#ifndef CVM_OPTIONS_HPP
#define CVM_OPTIONS_HPP
#include <memory>
#include "cvm/selection.hpp"
#include "cvm/enums.hpp"
namespace CVM
{
class Options
{
public:
    /// @brief Construtor
    Options(); 
    /// @brief Copy constructor.
    Options(const Options &options);
    /// @brief Move constructor.
    Options(Options &&options) noexcept;

    /// @brief Copy assignment operator.
    /// @param[in] options   The options class to copy to this.
    /// @result A deep copy of the input options.
    Options& operator=(const Options &options);
    /// @brief Move assignment operator.
    /// @param[in,out] options  The options class whose memory will be moved
    ///                         to this.  On exit, options's behavior is 
    ///                         undefined.
    /// @result The memory from options moved to this.
    Options& operator=(Options &&options) noexcept;

    /// @brief Parses and ini file.
    /// @param[in] fileName  The name of the ini file.
    void parse(const std::string &fileName);

    /// @brief Destructor
    ~Options();
    /// @brief Resets the class.
    void clear() noexcept;

    /// @name CVM
    /// @{
    /// @result The CVM model selection.
    Selection getSelection() const noexcept; 
    /// @result The P velocity file name for the given layer.
    /// @throws std::runtime_error if this was not read from the ini file.
    std::string getPVelocityFileName(LayerIdentifier layer) const;
    /// @result The S velocity file name for the given layer.
    /// @throws std::runtime_error if this was not read from the ini file.
    std::string getSVelocityFileName(LayerIdentifier layer) const;
    /// @result The velocity (m/s) to impute NaN P velocities with.
    ///         These values typically occur in water so by default this
    ///         is the speed of sound in water.
    double getPImputationVelocity() const;
    /// @result The velocity (m/s) to impute NaN S velocities with.  
    ///         These values typically occur in water so by default this is 0.
    double getSImputationVelocity() const;
    /// @}

    /// @name NonLinLoc
    /// @{
    /// @}
private:
    class OptionsImpl;
    std::unique_ptr<OptionsImpl> pImpl;
};
}
#endif
