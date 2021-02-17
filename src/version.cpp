#include <string>
#include "cvm/version.hpp"

using namespace CVM;

int Version::getMajor() noexcept
{
    return CVM_MAJOR;
}

int Version::getMinor() noexcept
{
    return CVM_MINOR;
}

int Version::getPatch() noexcept
{
    return CVM_PATCH;
}

bool Version::isAtLeast(const int major, const int minor,
                        const int patch) noexcept
{
    if (CVM_MAJOR < major){return false;}
    if (CVM_MAJOR > major){return true;}
    if (CVM_MINOR < minor){return false;}
    if (CVM_MINOR > minor){return true;}
    if (CVM_PATCH < patch){return false;}
    return true;
}

std::string Version::getVersion() noexcept
{
    std::string version(std::to_string(getMajor()) + "." 
                      + std::to_string(getMinor()) + "." 
                      + std::to_string(getPatch()));
    return version;
}
