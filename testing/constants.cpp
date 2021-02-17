#include <iostream>
#include <cmath>
#include <array>
#include <numeric>
#include "cvm/constants.hpp"
#include "cvm/geodetic.hpp"
#include <gtest/gtest.h>
using namespace CVM;
namespace
{
double getMaxLongitude();
double getMaxLatitude();
void getMaxNxNy(std::vector<int> *maxNx, std::vector<int> *maxNy);

TEST(Constants, Constants)
{
    Constants c;

    EXPECT_EQ(c.getNumberOfLayers(), 3);
    EXPECT_EQ(c.getUTMZone(), 10);

    EXPECT_EQ(c.getNumberOfGridPointsInX(LayerIdentifier::TOP),    3271);
    EXPECT_EQ(c.getNumberOfGridPointsInX(LayerIdentifier::MIDDLE), 2181); 
    EXPECT_EQ(c.getNumberOfGridPointsInX(LayerIdentifier::BOTTOM),  727);

    EXPECT_EQ(c.getNumberOfGridPointsInY(LayerIdentifier::TOP),    5367);
    EXPECT_EQ(c.getNumberOfGridPointsInY(LayerIdentifier::MIDDLE), 3578); 
    EXPECT_EQ(c.getNumberOfGridPointsInY(LayerIdentifier::BOTTOM), 1193);

    EXPECT_EQ(c.getNumberOfGridPointsInZ(LayerIdentifier::TOP),    13);
    EXPECT_EQ(c.getNumberOfGridPointsInZ(LayerIdentifier::MIDDLE), 29);
    EXPECT_EQ(c.getNumberOfGridPointsInZ(LayerIdentifier::BOTTOM), 55);

    EXPECT_NEAR(c.getGridSpacingInX(LayerIdentifier::TOP),    200, 1.e-13);
    EXPECT_NEAR(c.getGridSpacingInX(LayerIdentifier::MIDDLE), 300, 1.e-13);
    EXPECT_NEAR(c.getGridSpacingInX(LayerIdentifier::BOTTOM), 900, 1.e-13);

    EXPECT_NEAR(c.getGridSpacingInY(LayerIdentifier::TOP),    200, 1.e-13);
    EXPECT_NEAR(c.getGridSpacingInY(LayerIdentifier::MIDDLE), 300, 1.e-13);
    EXPECT_NEAR(c.getGridSpacingInY(LayerIdentifier::BOTTOM), 900, 1.e-13);

    EXPECT_NEAR(c.getGridSpacingInZ(LayerIdentifier::TOP),    100, 1.e-13);
    EXPECT_NEAR(c.getGridSpacingInZ(LayerIdentifier::MIDDLE), 300, 1.e-13);
    EXPECT_NEAR(c.getGridSpacingInZ(LayerIdentifier::BOTTOM), 900, 1.e-13);

    EXPECT_NEAR(c.getLayerStartDepth(LayerIdentifier::TOP),    0,     1.e-11);
    EXPECT_NEAR(c.getLayerStartDepth(LayerIdentifier::MIDDLE), 1500,  1.e-11);
    EXPECT_NEAR(c.getLayerStartDepth(LayerIdentifier::BOTTOM), 10800, 1.e-11);

    auto z0 = c.getLayerStartDepth(LayerIdentifier::TOP) 
            + (c.getNumberOfGridPointsInZ(LayerIdentifier::TOP) - 1)
             *c.getGridSpacingInZ(LayerIdentifier::TOP);
    auto z1 = c.getLayerStartDepth(LayerIdentifier::MIDDLE) 
            + (c.getNumberOfGridPointsInZ(LayerIdentifier::MIDDLE) - 1)
             *c.getGridSpacingInZ(LayerIdentifier::MIDDLE);
    auto z2 = c.getLayerStartDepth(LayerIdentifier::BOTTOM) 
            + (c.getNumberOfGridPointsInZ(LayerIdentifier::BOTTOM) - 1)
             *c.getGridSpacingInZ(LayerIdentifier::BOTTOM);
    EXPECT_NEAR(c.getLayerEndDepth(LayerIdentifier::TOP),    z0, 1.e-11);
    EXPECT_NEAR(c.getLayerEndDepth(LayerIdentifier::MIDDLE), z1, 1.e-11);
    EXPECT_NEAR(c.getLayerEndDepth(LayerIdentifier::BOTTOM), z2, 1.e-11);

    EXPECT_NEAR(c.getLayerStartDepth(LayerIdentifier::TOP),
                c.getMinimumDepth(), 1.e-12);
    EXPECT_NEAR(c.getLayerEndDepth(LayerIdentifier::BOTTOM),
                c.getMaximumDepth(), 1.e-12);

    auto mLCM = static_cast<double> (std::lcm(std::lcm(200, 300), 900));
    EXPECT_NEAR(c.getLeastCommonMultipleOfGridSpacingsInXAndY(), mLCM, 1.e-12);

    // Compute max usable lat/lon
    EXPECT_NEAR(c.getMaximumLatitude(),  getMaxLatitude(),  1.e-8);
    EXPECT_NEAR(c.getMaximumLongitude(), getMaxLongitude(), 1.e-8);

    //std::cout << std::setprecision(14) << c.getMaximumLatitude() << std::endl;
    //std::cout << std::setprecision(14) << c.getMaximumLongitude() << std::endl;
    std::vector<int> maxNx, maxNy;
    getMaxNxNy(&maxNx, &maxNy);
}

double getMaxLongitude()
{
    Constants c;
    auto lcm = c.getLeastCommonMultipleOfGridSpacingsInXAndY();
    double x1 = 0;
    for (int i = 0; i < c.getNumberOfLayers(); ++i)
    {
        auto layer = static_cast<LayerIdentifier> (i);
        auto nx = c.getNumberOfGridPointsInX(layer);
        auto dx = c.getGridSpacingInX(layer);
        auto xmax = static_cast<double> (nx - 1)*dx;
        auto nxLCM = static_cast<int> (xmax/lcm) + 1;
        x1 = std::max(x1, lcm*(nxLCM - 1)); 
    }
    double y1 = 0;
    for (int i = 0; i < c.getNumberOfLayers(); ++i)
    {
        auto layer = static_cast<LayerIdentifier> (i);
        auto ny = c.getNumberOfGridPointsInY(layer);
        auto dy = c.getGridSpacingInY(layer);
        auto ymax = static_cast<double> (ny - 1)*dy;
        auto nyLCM = static_cast<int> (ymax/lcm) + 1;
        y1 = std::max(y1, lcm*(nyLCM - 1));
    }
    // Longitude gets really messed up because UTM zone is wrong.
    // In this case, evaluate at the maximum latitude.
    x1 = x1 + c.getUTMOriginInX();
    y1 = y1 + c.getUTMOriginInY();
    //auto y0 = getUTMOriginInY();
    auto latLon = Geodetic::utmToLatitudeLongitude(std::pair(x1, y1),
                                                   c.getUTMZone());
//std::cout << std::setprecision(12) << "ll " << latLon.first << " " << latLon.second << std::endl;
    return latLon.second;
}

double getMaxLatitude()
{
    Constants c;
    auto lcm = c.getLeastCommonMultipleOfGridSpacingsInXAndY();
    double y1 = 0;
    for (int i = 0; i < c.getNumberOfLayers(); ++i)
    {
        auto layer = static_cast<LayerIdentifier> (i);
        auto ny = c.getNumberOfGridPointsInY(layer);
        auto dy = c.getGridSpacingInY(layer);
        auto ymax = static_cast<double> (ny - 1)*dy;
        auto nyLCM = static_cast<int> (ymax/lcm) + 1;
        y1 = std::max(y1, lcm*(nyLCM - 1));
    }
    auto x0 = c.getUTMOriginInX();
    y1 = y1 + c.getUTMOriginInY();
    auto latLon = Geodetic::utmToLatitudeLongitude(std::pair(x0, y1),
                                                   c.getUTMZone()); 
//std::cout << std::setprecision(12) << "llMaxLat " << latLon.first << " " << latLon.second << std::endl;
    return latLon.first;
}

void getMaxNxNy(std::vector<int> *maxNx, std::vector<int> *maxNy)
{
    Constants c;
    auto maxLat = c.getMaximumLatitude();
    auto maxLon = c.getMaximumLongitude();
    //auto maxLat = 49.770067370537;
    //auto maxLon = -121.04460751046;
    auto [utmX1, utmY1]
        = Geodetic::latitudeLongitudeToUTM(std::pair(maxLat, maxLon),
                                           c.getUTMZone()); 
    auto x0 = c.getUTMOriginInX();
    auto y0 = c.getUTMOriginInY();
    auto lcm = c.getLeastCommonMultipleOfGridSpacingsInXAndY();
    maxNx->resize(c.getNumberOfLayers(), 0);
    maxNy->resize(c.getNumberOfLayers(), 0);
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
                //if (ix + idx < nx){ixMax = ix + idx;}
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
                //if (iy + idy < ny){iyMax = iy + idy;}
                break;
            }
        }
        assert(ixMax > 0);
        assert(iyMax > 0);
        assert(ixMax <= c.getNumberOfGridPointsInX(layer));
        assert(iyMax <= c.getNumberOfGridPointsInY(layer));
        std::cout << ixMax << " " << iyMax << std::endl;
        maxNx->at(i) = ixMax;
        maxNy->at(i) = iyMax;
        auto xMax = x0 + dx*ixMax;
        auto yMax = y0 + dy*iyMax;
        auto [lat1, lon1]
            = Geodetic::utmToLatitudeLongitude(std::pair(xMax, yMax),
                                               c.getUTMZone());
        assert(lat1 <= c.getMaximumLatitude());
        assert(lon1 <= c.getMaximumLongitude());
        std::cout << std::setprecision(14) << lat1 << " " << lon1 << std::endl;
        std::cout << c.getMaximumLatitude() << " " << c.getMaximumLongitude() << std::endl;
    }
}

}
