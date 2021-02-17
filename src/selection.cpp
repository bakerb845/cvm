#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <cassert>
#include "cvm/selection.hpp"
#include "cvm/geodetic.hpp"
#include "cvm/constants.hpp"

using namespace CVM;

namespace
{
 
size_t layerToIndex(const LayerIdentifier &layer) noexcept
{
    return static_cast<size_t> (layer);
}

/*
int getMinimum(const int nx, const double x0, const double dx,
               const double xTarget, const int idx)
{
    int ix0 = 0;
    bool lfound = false;
    for (int ix = 0; ix < nx; ix = ix + idx)
    {   
        auto x = x0 + dx*static_cast<double> (ix);
        if (x > xTarget)
        {
            ix0 = std::max(0, ix - idx);
            lfound = true;
            break;
        }
     }
     if (!lfound){std::cerr << "getMinimum failed to bracket" << std::endl;}
#ifndef NDEBUG
     assert(ix0 >= 0 && ix0 < nx);
#endif
     return ix0;
}

int getMaximum(const int nx, const double x0, const double dx, 
               const double xTarget, const int idx)
{
    int ix1 = nx;
    bool lfound = false;
    for (int ix = 0; ix < nx; ix = ix + idx)
    {
        auto x = x0 + dx*static_cast<double> (ix);
        if (x >= xTarget)
        {
            ix1 = std::min(ix, nx);
            lfound = true;
            break;
        }
    }
    if (!lfound){std::cerr << "getMaximum failed to bracket" << std::endl;}
#ifndef NDEBUG
     assert(ix1 >= 0 && ix1 <= nx); 
#endif
     return ix1;
}
*/

void getMinNxNy(const double minLat, const double minLon,
                std::array<int, 3> *maxNx, std::array<int, 3> *maxNy)
{
    Constants c;
    auto [utmX1, utmY1]
        = Geodetic::latitudeLongitudeToUTM(std::pair(minLat, minLon),
                                           c.getUTMZone());
    auto x0 = c.getUTMOriginInX();
    auto y0 = c.getUTMOriginInY();
    auto lcm = c.getLeastCommonMultipleOfGridSpacingsInXAndY();
    //maxNx->resize(c.getNumberOfLayers(), 0); 
    //maxNy->resize(c.getNumberOfLayers(), 0); 
    for (int i = 0; i < c.getNumberOfLayers(); ++i)
    {
        auto layer = static_cast<LayerIdentifier> (i);
        auto nx = c.getNumberOfGridPointsInX(layer);
        auto ny = c.getNumberOfGridPointsInY(layer);
        auto dx = c.getGridSpacingInX(layer);
        auto dy = c.getGridSpacingInY(layer);
        auto idx = static_cast<int> (std::round(lcm/dx));
        auto idy = static_cast<int> (std::round(lcm/dy));
        int ixMax = 0;
        for (int ix = 0; ix < nx; ix = ix + idx)
        {
            auto x = x0 + dx*(ix + idx);
            if (x > utmX1)
            {
                ixMax = ix; 
                break;
            }   
        }
        int iyMax = 1;
        for (int iy = 0; iy < ny; iy = iy + idy)
        {
            auto y = y0 + dy*(iy + idy);
            if (y > utmY1)
            {
                iyMax = iy; 
                break;
            }   
        }   
#ifndef NDEBUG
        assert(ixMax > 0); 
        assert(iyMax > 0); 
        assert(ixMax <= c.getNumberOfGridPointsInX(layer));
        assert(iyMax <= c.getNumberOfGridPointsInY(layer));
#endif
        maxNx->at(i) = ixMax;
        maxNy->at(i) = iyMax;
#ifndef NDEBUG
        auto xMax = x0 + dx*ixMax;
        auto yMax = y0 + dy*iyMax;
        auto [lat1, lon1]
            = Geodetic::utmToLatitudeLongitude(std::pair(xMax, yMax),
                                               c.getUTMZone());
        assert(lat1 <= c.getMaximumLatitude());
        assert(lon1 <= c.getMaximumLongitude());
#endif
        //std::cout << std::setprecision(14) << lat1 << " " << lon1 << std::endl;
        //std::cout << c.getMaximumLatitude() << " " << c.getMaximumLongitude() << std::endl;
    }   
}


void getMaxNxNy(const double maxLat, const double maxLon,
                std::array<int, 3> *maxNx, std::array<int, 3> *maxNy)
{
    Constants c;
    auto [utmX1, utmY1]
        = Geodetic::latitudeLongitudeToUTM(std::pair(maxLat, maxLon),
                                           c.getUTMZone());
    auto x0 = c.getUTMOriginInX();
    auto y0 = c.getUTMOriginInY();
    auto lcm = c.getLeastCommonMultipleOfGridSpacingsInXAndY();
    //maxNx->resize(c.getNumberOfLayers(), 0);
    //maxNy->resize(c.getNumberOfLayers(), 0);
    for (int i = 0; i < c.getNumberOfLayers(); ++i)
    {
        auto layer = static_cast<LayerIdentifier> (i);
        auto nx = c.getNumberOfGridPointsInX(layer);
        auto ny = c.getNumberOfGridPointsInY(layer);
        auto dx = c.getGridSpacingInX(layer);
        auto dy = c.getGridSpacingInY(layer);
        auto idx = static_cast<int> (std::round(lcm/dx));
        auto idy = static_cast<int> (std::round(lcm/dy));
        int ixMax = 0;
        for (int ix = 0; ix < nx; ix = ix + idx)
        {
            auto x = x0 + dx*(ix + idx);
            if (x > utmX1)
            {
                ixMax = ix;
                if (ix + idx < c.getMaximumUsableNumberOfGridPointsInX(layer))
                {
                    ixMax = ix + idx;
                }
                break;
            }
        }
        int iyMax = 1;
        for (int iy = 0; iy < ny; iy = iy + idy)
        {
            auto y = y0 + dy*(iy + idy);
            if (y > utmY1)
            {
                iyMax = iy;
                if (iy + idy < c.getMaximumUsableNumberOfGridPointsInX(layer))
                {
                    iyMax = iy + idy;
                }
                break;
            }
        }
#ifndef NDEBUG
        assert(ixMax > 0);
        assert(iyMax > 0);
        assert(ixMax <= c.getNumberOfGridPointsInX(layer));
        assert(iyMax <= c.getNumberOfGridPointsInY(layer));
#endif
        maxNx->at(i) = ixMax;
        maxNy->at(i) = iyMax;
#ifndef NDEBUG
        auto xMax = x0 + dx*ixMax;
        auto yMax = y0 + dy*iyMax;
        auto [lat1, lon1]
            = Geodetic::utmToLatitudeLongitude(std::pair(xMax, yMax),
                                               c.getUTMZone());
        assert(lat1 <= c.getMaximumLatitude());
        assert(lon1 <= c.getMaximumLongitude());
std::cout << lat1 << " " << lon1 << std::endl;
#endif
        //std::cout << std::setprecision(14) << lat1 << " " << lon1 << std::endl;
        //std::cout << c.getMaximumLatitude() << " " << c.getMaximumLongitude() << std::endl;
    }
}



}

class Selection::SelectionImpl
{
public:
    void updateGridIndices()
    {
        if (mReadAll)
        {
            mEndX = std::array<int, 3>{
               Constants::getMaximumUsableNumberOfGridPointsInX(LayerIdentifier::TOP),
               Constants::getMaximumUsableNumberOfGridPointsInX(LayerIdentifier::MIDDLE),
               Constants::getMaximumUsableNumberOfGridPointsInX(LayerIdentifier::BOTTOM)
            };
            mStartX = std::array<int, 3>{0, 0, 0};
            mEndY = std::array<int, 3>{
               Constants::getMaximumUsableNumberOfGridPointsInY(LayerIdentifier::TOP),
               Constants::getMaximumUsableNumberOfGridPointsInY(LayerIdentifier::MIDDLE),
               Constants::getMaximumUsableNumberOfGridPointsInY(LayerIdentifier::BOTTOM)
            };
            mStartZ = std::array<int, 3>{0, 0, 0};
            mEndZ = {
               Constants::getNumberOfGridPointsInZ(LayerIdentifier::TOP),
               Constants::getNumberOfGridPointsInZ(LayerIdentifier::MIDDLE),
               Constants::getNumberOfGridPointsInZ(LayerIdentifier::BOTTOM)
            };
            mUpdateGridPoints = false;
            return;
        }
std::cout << mMinimumLatitude << " " << mMaximumLatitude << std::endl;
std::cout << mMinimumLongitude << " " << mMaximumLongitude << std::endl;
        getMinNxNy(mMinimumLatitude, mMinimumLongitude, &mStartX, &mStartY);
        getMaxNxNy(mMaximumLatitude, mMaximumLongitude, &mEndX,   &mEndY);
        mUpdateGridPoints = false;
        /*
        auto [utmx0, utmy0] = Geodetic::latitudeLongitudeToUTM(
            std::pair(mMinimumLatitude, mMinimumLongitude),
            Constants::getUTMZone());
        auto [utmx1, utmy1] = Geodetic::latitudeLongitudeToUTM(
            std::pair(mMaximumLatitude, mMaximumLongitude),
            Constants::getUTMZone());
        auto x0 = Constants::getUTMOriginInX();
        auto y0 = Constants::getUTMOriginInY();
        auto d = Constants::getLeastCommonMultipleOfGridSpacingsInXAndY();
        for (int layer = 0; layer < Constants::getNumberOfLayers(); ++layer)
        {
            auto layerID = static_cast<LayerIdentifier> (layer);
            auto dx = Constants::getGridSpacingInX(layerID);
            auto dy = Constants::getGridSpacingInY(layerID);
            auto incDx = static_cast<int> (std::round(d/dx));
            auto incDy = static_cast<int> (std::round(d/dy));
            auto nx = Constants::getMaximumUsableNumberOfGridPointsInX(layerID);
            auto ny = Constants::getMaximumUsableNumberOfGridPointsInY(layerID);
            // Get minimum and maximum indices
std::cout << x0 << " " << y0 << " " << utmx0 << " " << utmy0 << " "
          << mMinimumLatitude << " " << mMinimumLongitude << " " << mMaximumLatitude << " " << mMaximumLongitude << std::endl;
            mStartX[layer] = getMinimum(nx, x0, dx, utmx0, incDx);
            mStartY[layer] = getMinimum(ny, y0, dy, utmy0, incDy);
            mEndX[layer] = getMaximum(nx, x0, dx, utmx1, incDx);
            mEndY[layer] = getMaximum(ny, y0, dy, utmy1, incDy);
#ifndef NDEBUG
            assert(mStartX[layer] < mEndX[layer]);
            assert(mStartY[layer] < mEndY[layer]);
#endif 
        }
        */
        mReadAll = false;
    }
//private:
    std::array<int, 3> mStartX{0, 0, 0};
    std::array<int, 3> mEndX{
      Constants::getMaximumUsableNumberOfGridPointsInX(LayerIdentifier::TOP),
      Constants::getMaximumUsableNumberOfGridPointsInX(LayerIdentifier::MIDDLE),
      Constants::getMaximumUsableNumberOfGridPointsInX(LayerIdentifier::BOTTOM)
    };
    std::array<int, 3> mStartY{0, 0, 0};
    std::array<int, 3> mEndY{
      Constants::getMaximumUsableNumberOfGridPointsInY(LayerIdentifier::TOP),
      Constants::getMaximumUsableNumberOfGridPointsInY(LayerIdentifier::MIDDLE),
      Constants::getMaximumUsableNumberOfGridPointsInY(LayerIdentifier::BOTTOM)
    }; 
    std::array<int, 3> mStartZ{0, 0, 0};
    std::array<int, 3> mEndZ{
        Constants::getNumberOfGridPointsInZ(LayerIdentifier::TOP),
        Constants::getNumberOfGridPointsInZ(LayerIdentifier::MIDDLE),
        Constants::getNumberOfGridPointsInZ(LayerIdentifier::BOTTOM)}; 
    double mMinimumLatitude  = Constants::getMinimumLatitude();
    double mMinimumLongitude = Constants::getMinimumLongitude(); 
    double mMaximumLatitude  = Constants::getMaximumLatitude();
    double mMaximumLongitude = Constants::getMaximumLongitude();
    bool mReadAll = true;
    bool mUpdateGridPoints = true;
};

/// C'tor
Selection::Selection() :
    pImpl(std::make_unique<SelectionImpl> ())
{
}

/// Copy c'tor
Selection::Selection(const Selection &selection)
{
    *this = selection;
}

/// Move c'tor
Selection::Selection(Selection &&selection) noexcept
{
    *this = std::move(selection);
}

/// Copy assignment
Selection& Selection::operator=(const Selection &selection)
{
    if (&selection == this){return *this;}
    pImpl = std::make_unique<SelectionImpl> (*selection.pImpl);
    return *this;
}

/// Move assignment
Selection& Selection::operator=(Selection &&selection) noexcept
{
    if (&selection == this){return *this;}
    pImpl = std::move(selection.pImpl);
    return *this;
}

/// Destructor
Selection::~Selection() = default;

/// Load entire model?
bool Selection::readAll() const noexcept
{
    return pImpl->mReadAll;
}

/// Set minimum latitude and longitude
void Selection::setMinimumLatitudeAndLongitude(
    const std::pair<double, double> latLonIn)
{
    auto latlon = latLonIn;
    if (latlon.second > 180)
    {
        latlon.second = latlon.second - 360;
    }
    if (latlon.first < Constants::getMinimumLatitude())
    {
        throw std::invalid_argument("latitude must be at least "
                             + std::to_string(Constants::getMinimumLatitude()));
    }
    if (latlon.second < Constants::getMinimumLongitude())
    {   
        throw std::invalid_argument("longitude must be at least "
                            + std::to_string(Constants::getMinimumLongitude()));
    }
    if (latlon.first > pImpl->mMaximumLatitude)
    {
        std::cout << "Resetting maximum latitude to CVM maximum latitude"
                  << std::endl;
        pImpl->mMaximumLatitude = Constants::getMaximumLatitude();
    }
    if (latlon.second > pImpl->mMaximumLongitude)
    {
        std::cout << "Resetting maximum longitude to CVM maximum longitude"
                  << std::endl;
        pImpl->mMaximumLongitude = Constants::getMaximumLongitude();
    }
    pImpl->mMinimumLatitude = latlon.first;
    pImpl->mMinimumLongitude = latlon.second;
    pImpl->mUpdateGridPoints = true;
    pImpl->mReadAll = false;
    //pImpl->updateGridIndices();
}

/// Set maximum latitude and longitude
void Selection::setMaximumLatitudeAndLongitude(
    const std::pair<double, double> latLonIn)
{
    auto latlon = latLonIn;
    if (latlon.second > 180)
    {
         latlon.second = latlon.second - 360;
    }
    if (latlon.first > Constants::getMaximumLatitude())
    {
        throw std::invalid_argument("latitude cannot exceed "
                              + std::to_string(Constants::getMaximumLatitude()));
    }
    if (latlon.second > Constants::getMaximumLongitude())
    {   
        throw std::invalid_argument("longitude cannot exceed "
                            + std::to_string(Constants::getMaximumLongitude()));
    }
    if (latlon.first < pImpl->mMinimumLatitude)
    {
        std::cout << "Resetting minimum latitude to model origin" << std::endl;
        pImpl->mMinimumLatitude = Constants::getMinimumLatitude();
    }
    if (latlon.second < pImpl->mMinimumLongitude)
    {
        std::cout << "Resetting minimum longitude to model origin" << std::endl;
        pImpl->mMinimumLongitude = Constants::getMinimumLongitude();
    }
    pImpl->mMaximumLatitude = latlon.first;
    pImpl->mMaximumLongitude = latlon.second;
    pImpl->mUpdateGridPoints = true;
    pImpl->mReadAll = false;
    //pImpl->updateGridIndices();
}

/// Get start x index
int Selection::getStartPointInX(LayerIdentifier layer) const noexcept
{
    if (pImpl->mUpdateGridPoints)
    {
        pImpl->updateGridIndices();
    }
    return pImpl->mStartX[layerToIndex(layer)];
}

/// Get end x index
int Selection::getEndPointInX(LayerIdentifier layer) const noexcept
{
    if (pImpl->mUpdateGridPoints)
    {
        pImpl->updateGridIndices();
    }
    return pImpl->mEndX[layerToIndex(layer)];
}

/// Get start y index
int Selection::getStartPointInY(LayerIdentifier layer) const noexcept
{
    if (pImpl->mUpdateGridPoints)
    {   
        pImpl->updateGridIndices();
    }
    return pImpl->mStartY[layerToIndex(layer)];
}
 
/// Get end y index
int Selection::getEndPointInY(LayerIdentifier layer) const noexcept
{
    if (pImpl->mUpdateGridPoints)
    {   
        pImpl->updateGridIndices();
    }
    return pImpl->mEndY[layerToIndex(layer)];
}

