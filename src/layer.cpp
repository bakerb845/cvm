#include <iostream>
#include <bit>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#ifndef NDEBUG
  #include <cassert>
#endif
#ifdef USE_STD_FILESYSTEM
  #include <filesystem>
  namespace fs = std::filesystem;
  #define HAVE_FS 1
#endif
#ifdef USE_STD_EXPERIMENTAL_FILESYSTEM
  #include <experimental/filesystem>
  namespace = fs = std::experimental::filesystem;
  #define HAVE_FS 1
#endif
#include "cvm/layer.hpp"
#include "cvm/geodetic.hpp"
#include "cvm/selection.hpp"
#include "cvm/options.hpp"
#include "cvm/constants.hpp"
#include "private/interpolate.hpp"

#define VP_WATER 1482.0
#define VS_WATER 0.0
#define NAN_CVM 1.e20

using namespace CVM;

namespace
{

/// @brief Converts from the Vp to the density.
template<typename T> 
T densityBrocher(T Vp)
{
    T rho, Vp2, Vp3, Vp4, Vp5;
    const T c1 = 1.6612;
    const T c2 =-0.4721;
    const T c3 = 0.0671;
    const T c4 =-0.0043;
    const T c5 = 0.000106;
    Vp2 = Vp*Vp;
    Vp3 = Vp*Vp2;
    Vp4 = Vp*Vp3;
    Vp5 = Vp*Vp4; 
    rho = c1*Vp + c2*Vp2 + c3*Vp3 + c4*Vp4 + c5*Vp5;
    return rho;
}

/*
/// Reads a binary file
std::vector<char> readBinaryFile(const std::string &fileName)
{
    // Read the binary file
    std::ifstream cvmFile(fileName, std::ios::in | std::ios::binary);
    cvmFile.seekg(0, cvmFile.end);
    size_t nbytes = cvmFile.tellg();  
    if (nbytes == 0)
    {
        throw std::invalid_argument(fileName + " has no data");
    }
    // Return to start of file 
    std::vector<char> result;
    cvmFile.seekg(0, cvmFile.beg);
    // Now read it
    std::vector<char> buffer(nbytes);
    cvmFile.read(buffer.data(), buffer.size()*sizeof(char));
    return buffer; 
}
*/

}

template<LayerIdentifier E>
class Layer<E>::LayerImpl
{
public:
    // Get a pointer to part of the data
    float *getVelocityPointer(int ix, int iy, int iz, const bool isP)
    {
#ifndef NDEBUG
        assert(ix >= 0 && ix < mNx);
        assert(iy >= 0 && iy < mNy);
        assert(iz >= 0 && iz < mNz);
#endif
        auto indx = static_cast<size_t> (iz*mNy*mNx + iy*mNx + ix);
        if (isP)
        {
            return mVp.data() + indx;
        }
        else
        {
            return mVs.data() + indx;
        }
    }
    // Clear memory for a Vp or Vs vector
    void clearVelocity(const bool isP)
    {
        if (isP)
        {
            mVp.clear();
            mRho.clear();
            mHavePVelocity = false;
            mHaveDensity = false;
        }
        else
        {
            mVs.clear();
            mHaveSVelocity = false;
        }
    }
    // Resize Vp or Vs vector
    void resizeVelocity(const size_t nxyz, const bool isP)
    {
        if (isP)
        {
            mVp.resize(nxyz, 0);
            mHavePVelocity = false;
            mHaveDensity = false;
        }
        else
        {
            mVs.resize(nxyz, 0);
            mHaveSVelocity = false;
        }
    }
    std::vector<float> mVp;
    std::vector<float> mVs;
    std::vector<float> mRho;
    std::pair<double, double> mLowerLeftUTM;
    std::pair<double, double> mUpperRightUTM;
    int mNx = 0;
    int mNy = 0;
    int mNz = 0;
    bool mHavePVelocity = false;
    bool mHaveSVelocity = false;
    bool mHaveDensity = false;
};

/// Release memory and reset the class
template<LayerIdentifier E>
void Layer<E>::clear() noexcept
{
    pImpl->mVp.clear();
    pImpl->mVs.clear();
    pImpl->mRho.clear();
    pImpl->mHavePVelocity = false;
    pImpl->mHaveSVelocity = false;
    pImpl->mHaveDensity = false;
}

/// C'tor
template<LayerIdentifier E>
Layer<E>::Layer() :
    pImpl(std::make_unique<LayerImpl> ())
{
}

/// Destructor
template<LayerIdentifier E>
Layer<E>::~Layer() = default;

/// Load the layer
template<LayerIdentifier E>
void Layer<E>::load(const Options &options, const bool isP)
{
    std::string fileName;
    if (isP)
    {
        std::cout << "Loading Vp in layer: "
                  << std::to_string(static_cast<int> (E) + 1) << std::endl;
        fileName = options.getPVelocityFileName(E);
    }
    else
    {
        std::cout << "Loading Vs in layer: "
                  << std::to_string(static_cast<int> (E) + 1) << std::endl;
        fileName = options.getSVelocityFileName(E);
    }
    auto selection = options.getSelection();
#ifdef HAVE_FS
    if (!fs::exists(fileName))
    {
        throw std::invalid_argument("File: " + fileName + " does not exist");
    }
#endif
    // Reset velocity in layer
    pImpl->clearVelocity(isP);
    // Get nx, ny, nz for this layer
    Constants constants;
    auto nx = constants.getNumberOfGridPointsInX(E);
    auto ny = constants.getNumberOfGridPointsInY(E);
    auto nz = constants.getNumberOfGridPointsInZ(E);
    auto nBytesExpected = static_cast<size_t> (nx*ny*nz)*sizeof(float);
    // Am I swapping bytes?
    bool lSwap = false;
    if (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
    {
        lSwap = true;
        std::cout << "Will swap bytes - warning this isn't tested" << std::endl;
    }
    else if (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
    {
        lSwap = false;
    }
    else
    {
        throw std::runtime_error("Don't know how to handle mixed endian");
    }
    // Ensure the binary file size makes sense
    std::ifstream cvmFile(fileName, std::ios::in | std::ios::binary);
    cvmFile.seekg(0, cvmFile.end);
    size_t nBytes = cvmFile.tellg();
    if (nBytes != nBytesExpected)
    {   
        if (nBytes == 0)
        {
            throw std::invalid_argument(fileName + " has no data");
        }
        throw std::invalid_argument(fileName + " has unexpected size");
    }
    // Go back to beginning of file
    cvmFile.seekg(0, cvmFile.beg);
    bool readEntireFile = selection.readAll();
    // Load entire file
    if (readEntireFile)
    {
        pImpl->mNx = nx;
        pImpl->mNy = ny;
        pImpl->mNz = nz;
        auto nxyzRead = static_cast<size_t> (nx*ny*nz);
        pImpl->resizeVelocity(nxyzRead, isP);
        if (!lSwap)
        {
            // Simply need a pointer to preallocated space then read it
            auto cBufferPtr = reinterpret_cast<char *>
                              (pImpl->getVelocityPointer(0, 0, 0, isP));
            cvmFile.read(cBufferPtr, nBytes);
        } 
        else
        {
            // Allocate space and load
            std::vector<char> cBuffer(nBytes);
            cvmFile.read(cBuffer.data(), nBytes);
            // Swap bytes
            for (size_t i = 0; i < nxyzRead; ++i)
            {
                std::reverse(cBuffer.data() + 4*i,  cBuffer.data() + 4*(i + 1));
            }
            // Copy result
            auto vPtrRead = reinterpret_cast<float *> (cBuffer.data()); 
            auto vPtr = pImpl->getVelocityPointer(0, 0, 0, isP);
            std::copy(vPtrRead, vPtrRead + nxyzRead, vPtr);
        }
    }
    else // Load a selection
    {
        auto utmX0 = Constants::getUTMOriginInX()
                   + Constants::getGridSpacingInX(E)
                    *selection.getStartPointInX(E);

        auto utmY0 = Constants::getUTMOriginInY()
                   + Constants::getGridSpacingInY(E)
                    *selection.getStartPointInY(E);

        auto utmX1 = Constants::getUTMOriginInX()
                   + Constants::getGridSpacingInX(E)
                    *(selection.getEndPointInX(E));
        auto utmY1 = Constants::getUTMOriginInY()
                   + Constants::getGridSpacingInY(E)
                    *(selection.getEndPointInY(E));
        auto [lat0, lon0] = Geodetic::utmToLatitudeLongitude(std::pair(utmX0, utmY0));
        auto [lat1, lon1] = Geodetic::utmToLatitudeLongitude(std::pair(utmX1, utmY1));
        std::cout << "Extracting: (" << lat0 << "," << lon0 << ") to (" 
                  << lat1 << "," << lon1 << ")" << " in layer " 
                  << static_cast<int> (E) + 1 << std::endl;
        pImpl->mLowerLeftUTM  = std::pair(utmX0, utmY0);
        pImpl->mUpperRightUTM = std::pair(utmX1, utmY1);          

        int ix0 = selection.getStartPointInX(E); //0;
        int ix1 = selection.getEndPointInX(E) + 1; //1;
        int iy0 = selection.getStartPointInY(E); //0;
        int iy1 = selection.getEndPointInY(E) + 1; //1;
        int iz0 = 0;
        int iz1 = nz;
        // Update
        pImpl->mNx = ix1 - ix0;
        pImpl->mNy = iy1 - iy0;
        pImpl->mNz = iz1 - iz0;
        auto nxyzRead
            = static_cast<size_t> ((ix1 - ix0)*(iy1 - iy0)*(iz1 - iz0));
        pImpl->resizeVelocity(nxyzRead, isP);
        int nxRead = ix1 - ix0 + 1;
        auto nBytesRead = static_cast<size_t> (nxRead)*sizeof(float);
        char *cBufferPtr = nullptr;
        std::vector<char> cBuffer;
        if (lSwap)
        {
            cBuffer.resize(nBytesRead);
            cBufferPtr = cBuffer.data();
        }
        // Load the file
        for (int iz = iz0; iz < iz1; ++iz)
        {
            for (int iy = iy0; iy < iy1; ++iy)
            {
                // Update pointer if reading directly to internal memory
                if (!lSwap)
                {
                    cBufferPtr
                        = reinterpret_cast<char *>
                          (pImpl->getVelocityPointer(0, iy-iy0, iz-iz0, isP));
                }
                // Find position in CVM.  Have to flip z.
                auto jz = nz - 1 - iz; 
                auto i1 = static_cast<size_t> (nx*ny*jz + nx*iy + ix0)
                         *sizeof(float);
                cvmFile.seekg(i1, cvmFile.beg);
                cvmFile.read(cBufferPtr, nBytesRead);
                // Swap bytes if necessary
                if (lSwap)
                {
                    for (int ix = 0; ix < pImpl->mNx; ++ix)
                    {
                        std::reverse(cBufferPtr + ix*4,
                                     cBufferPtr + (ix + 1)*4);
                    }  
                    auto vPtr = pImpl->getVelocityPointer(0, iy-iy0, iz-iz0,
                                                          isP); 
                    auto vPtrRead = reinterpret_cast<float *> (cBufferPtr);
                    std::copy(vPtrRead, vPtrRead + pImpl->mNx, vPtr);
                }
            } // Loop on y
            //std::cout << pImpl->mVp[iz] << std::endl;
        } // Loop on z
    }
    // Done with file
    cvmFile.close();
    // Fix water problem
    if (isP)
    {
        auto vpWater = static_cast<float> (options.getPImputationVelocity());
        std::replace_if(pImpl->mVp.begin(), pImpl->mVp.end(),
                        [&](auto x){return x > NAN_CVM;},
                        vpWater); 
        auto [vpMin, vpMax] = std::minmax_element(pImpl->mVp.begin(),
                                                  pImpl->mVp.end());
        std::cout << "Min/max imputed Vp: "
                  << *vpMin << " " << *vpMax << std::endl;
    }
    else
    {
        auto vsWater = static_cast<float> (options.getSImputationVelocity());
        std::replace_if(pImpl->mVs.begin(), pImpl->mVs.end(),
                        [&](auto x){return x > NAN_CVM;},
                        vsWater); 
        auto [vsMin, vsMax] = std::minmax_element(pImpl->mVs.begin(),
                                                  pImpl->mVs.end());
        std::cout << "Min/max imputed Vs: "
                  << *vsMin << " " << *vsMax << std::endl;
    }
    //for (auto vp : pImpl->mVp)
    //{
    //    std::cout << vp << " "; // << std::endl;
    //}
    // It's loaded
    if (isP)
    {
        pImpl->mHavePVelocity = true;
    }
    else
    {
        pImpl->mHaveSVelocity = true;
    }

    std::cout << std::endl;
}

/// Gets the layer identifier
template<LayerIdentifier E>
LayerIdentifier Layer<E>::getLayerIdentifier() const noexcept
{
    return E;
}

/// Get a pointer to the Vp
template<LayerIdentifier E>
const float *Layer<E>::getPVelocityPointer() const
{
    if (!havePVelocity())
    {
        throw std::runtime_error("P velocity in layer not yet loaded");
    }
    return pImpl->mVp.data();
} 

/// Get a pointer to the Vs
template<LayerIdentifier E>
const float *Layer<E>::getSVelocityPointer() const
{
    if (!haveSVelocity())
    {   
        throw std::runtime_error("S velocity in layer not yet loaded");
    }   
    return pImpl->mVs.data();
} 

/// Get a pointer to the density 
template<LayerIdentifier E>
const float *Layer<E>::getDensityPointer() const
{
    if (!havePVelocity())
    {
        throw std::runtime_error("Layer not yet loaded");
    }   
    // Compute density if not already done so
    if (!pImpl->mHaveDensity)
    {
        auto n = pImpl->mVp.size();
        const float *__restrict__ vPtr = getPVelocityPointer();
        pImpl->mRho.resize(n);
        float *dPtr = pImpl->mRho.data();
        for (size_t i = 0; i < n; ++i)
        {
            dPtr[i] = densityBrocher(vPtr[i]);
        }
        pImpl->mHaveDensity = true;
    }
    return pImpl->mRho.data();
} 

/// Have I read the velocity information yet?
template<LayerIdentifier E>
bool Layer<E>::havePVelocity() const noexcept
{
    return pImpl->mHavePVelocity;
}

template<LayerIdentifier E>
bool Layer<E>::haveSVelocity() const noexcept
{
    return pImpl->mHaveSVelocity;
}

/// Model lower left/upper right
template<LayerIdentifier E>
std::pair<double, double> Layer<E>::getSouthWestCornerUTM() const
{
    if (!havePVelocity() && !haveSVelocity())
    {
        throw std::runtime_error("No velocity model read");
    }
    return pImpl->mLowerLeftUTM;
}

template<LayerIdentifier E>
std::pair<double, double> Layer<E>::getNorthEastCornerUTM() const
{
    if (!havePVelocity() && !haveSVelocity())
    {   
        throw std::runtime_error("No velocity model read");
    }   
    return pImpl->mUpperRightUTM;
}

/// Nx
template<LayerIdentifier E>
int Layer<E>::getNumberOfGridPointsInX() const
{
    if (!havePVelocity() && !haveSVelocity())
    {   
        throw std::runtime_error("No velocity model read");
    }
    return pImpl->mNx;
}

/// Ny
template<LayerIdentifier E>
int Layer<E>::getNumberOfGridPointsInY() const
{
    if (!havePVelocity() && !haveSVelocity())
    {
        throw std::runtime_error("No velocity model read");
    }
    return pImpl->mNy;
}

/// Nz
template<LayerIdentifier E>
int Layer<E>::getNumberOfGridPointsInZ() const
{
    if (!havePVelocity() && !haveSVelocity())
    {
        throw std::runtime_error("No velocity model read");
    }
    return pImpl->mNz;
}



///--------------------------------------------------------------------------///
///                          Template Instantiation                          ///
///--------------------------------------------------------------------------///
template class CVM::Layer<LayerIdentifier::TOP>;
template class CVM::Layer<LayerIdentifier::MIDDLE>;
template class CVM::Layer<LayerIdentifier::BOTTOM>;
