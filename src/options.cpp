#include <iostream>
#include <array>
#include <string>
#ifdef USE_STD_FILESYSTEM
  #include <filesystem>
#endif
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include "cvm/options.hpp"
#include "cvm/constants.hpp"
#include "cvm/selection.hpp"

using namespace CVM;

class Options::OptionsImpl
{
public:
    Selection mSelection;
    std::array<std::string, 3> mPModelName{"", "", ""};
    std::array<std::string, 3> mSModelName{"", "", ""};
};

/// C'tor
Options::Options() :
    pImpl(std::make_unique<OptionsImpl> ())
{
    if (pImpl->mPModelName.size() !=
        static_cast<size_t> (Constants::getNumberOfLayers()))
    {
        throw std::runtime_error(
            "Internal error - update options to match cvm constants");
    }
}

Options::Options(const Options &options)
{
    *this = options;
}

Options::Options(Options &&options) noexcept
{
    *this = std::move(options);
}

/// Assignment operators
Options& Options::operator=(const Options &options)
{
    if (&options == this){return *this;}
    pImpl = std::make_unique<OptionsImpl> (*options.pImpl);
    return *this;
}

Options& Options::operator=(Options &&options) noexcept
{
    if (&options == this){return *this;}
    pImpl = std::move(options.pImpl);
    return *this;
}

/// Destructor
Options::~Options() = default;

/// Reset the class
void Options::clear() noexcept
{
    for (size_t i = 0; i < pImpl->mPModelName.size(); ++i)
    {
        pImpl->mPModelName[i].clear();
        pImpl->mSModelName[i].clear();
    }
}

/// Parses an ini file
void Options::parse(const std::string &fileName)
{
    boost::property_tree::ptree pt;
    boost::property_tree::read_ini(fileName, pt);
    pImpl->mPModelName[0] = pt.get<std::string>("CVM.layer1PModelName");
    pImpl->mPModelName[1] = pt.get<std::string>("CVM.layer2PModelName"); 
    pImpl->mPModelName[2] = pt.get<std::string>("CVM.layer3PModelName");

    pImpl->mSModelName[0] = pt.get<std::string>("CVM.layer1SModelName");
    pImpl->mSModelName[1] = pt.get<std::string>("CVM.layer2SModelName"); 
    pImpl->mSModelName[2] = pt.get<std::string>("CVM.layer3SModelName");
#ifdef USE_STD_FILESYSTEM
    for (size_t i = 0; i < pImpl->mPModelName.size(); ++i)
    {
        if (!std::filesystem::exists(pImpl->mPModelName[i]))
        {
            std::cerr << "P model name: " << pImpl->mPModelName[i]
                      << " does not exist" << std::endl;
            pImpl->mPModelName[i].clear();
        }
        if (!std::filesystem::exists(pImpl->mSModelName[i]))
        {
            std::cerr << "S model name: " << pImpl->mSModelName[i]
                      << " does not exist" << std::endl;
            pImpl->mSModelName[i].clear();
        }
    }
#endif

    Selection selection;
    double lat0 = Constants::getMinimumLatitude(); 
    double lat1 = Constants::getMaximumLatitude();
    double lon0 = Constants::getMinimumLongitude(); 
    double lon1 = Constants::getMaximumLongitude();
    try
    {
        lat0 = pt.get<double> ("Selection.latitude0");
        lat1 = pt.get<double> ("Selection.latitude1");
        lon0 = pt.get<double> ("Selection.longitude0");
        lon1 = pt.get<double> ("Selection.longitude1");
        if (lon0 > 180){lon0 = lon0 - 360;}
        if (lon1 > 180){lon1 = lon1 - 360;}
    }
    catch (const std::exception &e)
    {
    }
    // Verify bounds on selection
    if (lat0 < Constants::getMinimumLatitude())
    {
        lat0 = Constants::getMinimumLatitude();
        std::cerr << "Overriding latitude0 to " << lat0 << std::endl;
    }
    if (lat1 > Constants::getMaximumLatitude())
    {
        lat1 = Constants::getMaximumLatitude();
        std::cerr << "Overriding latitude1 to " << lat1 << std::endl;
    }
    if (lon0 < Constants::getMinimumLongitude())
    {
        lon0 = Constants::getMinimumLongitude();
        std::cerr << "Overriding longitude0 to " << lon0 << std::endl;
    }
    if (lon1 > Constants::getMaximumLongitude())
    {
        lon1 = Constants::getMaximumLongitude();
        std::cerr << "Overriding longitude1 to " << lon1 << std::endl;
    }
    // Selection must satisfy some strict criteria
    if (lat1 < lat0)
    {
        throw std::invalid_argument("latitude1 = " + std::to_string(lat1)
                                  + " cannot be less than latitude0 = "
                                  + std::to_string(lat0));
    } 
    if (lon1 < lon0)
    {
        throw std::invalid_argument("longitude1 = " + std::to_string(lon1)
                                  + " cannot be less than longitude0 = "
                                  + std::to_string(lon0));
    }
    // Set the selection
    //std::cout << lat0 << " " << lon0 << std::endl;
    //std::cout << lat1 << " " << lon1 << std::endl;
    selection.setMinimumLatitudeAndLongitude(std::pair(lat0, lon0)); 
    selection.setMaximumLatitudeAndLongitude(std::pair(lat1, lon1));
    pImpl->mSelection = selection;

    double depth0 = Constants::getMinimumDepth();
    double depth1 = Constants::getMaximumDepth();
    try
    {
        double depth0 = pt.get<double> ("Selection.depth0");
        double depth1 = pt.get<double> ("Selection.depth1");
        depth0 = depth0/1000;
        depth1 = depth1/1000;
    }
    catch (const std::exception &e)
    {
    }
    if (depth1 < depth0)
    {
        throw std::invalid_argument("Max cannot be less than minimum depth");
    }
    if (depth0 < Constants::getMinimumDepth())
    {
        auto errmsg = "Minimum model depth must be at least "
                    + std::to_string(Constants::getMinimumDepth()/1000);
        throw std::invalid_argument(errmsg);
    }
    if (depth1 > Constants::getMaximumDepth())
    {
        auto errmsg = "Maximum model depth cannot exceed "
                    + std::to_string(Constants::getMaximumDepth()/1000);
        throw std::invalid_argument(errmsg); 
    }
    selection.setMinimumAndMaximumDepth(std::pair(depth0, depth1));
}

Selection Options::getSelection() const noexcept
{
    return pImpl->mSelection;
}

/// Get P velocity file name
std::string Options::getPVelocityFileName(const LayerIdentifier layer) const
{
    auto i = static_cast<size_t> (layer);
    if (pImpl->mPModelName[i].empty())
    {
        throw std::runtime_error("P file name not read from ini file");
    }
    return pImpl->mPModelName[i];
}

/// Get S velocity file name
std::string Options::getSVelocityFileName(const LayerIdentifier layer) const
{
    auto i = static_cast<size_t> (layer);
    if (pImpl->mSModelName[i].empty())
    {   
        throw std::runtime_error("S file name not read from ini file");
    }
    return pImpl->mSModelName[i];
}


