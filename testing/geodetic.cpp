#include <iostream>
#include "cvm/geodetic.hpp"
#include <gtest/gtest.h>
namespace
{
TEST(CVM, Geodetic)
{
    const double tol = 10; // 10 meters is more than close enough
    std::pair<double, double> seattle{47.6062, -122.3321};
    std::pair<double, double> seattleUTM{550200, 5272748.59};
    auto xy = CVM::Geodetic::latitudeLongitudeToUTM(seattle, 10);
    EXPECT_NEAR(xy.first,  seattleUTM.first,  tol);
    EXPECT_NEAR(xy.second, seattleUTM.second, tol);
    EXPECT_NEAR(xy.first,  550200.21335878095, 1.e-8); // From fortran code
    EXPECT_NEAR(xy.second, 5272751.6364032440, 1.e-7);
    auto ll = CVM::Geodetic::utmToLatitudeLongitude(xy, 10);
    EXPECT_NEAR(ll.first,  seattle.first,  1.e-4);
    EXPECT_NEAR(ll.second, seattle.second, 1.e-4);
    EXPECT_NEAR(ll.first,   47.606227487491026, 1.e-10); // From Fortran code
    EXPECT_NEAR(ll.second, -122.33209965126906, 1.e-10);

// std::cout << std::setprecision(14) << xy.first << " " << xy.second - seattleUTM.second << std::endl;
//std::cout << ll.first << " " << ll.second << std::endl;
//    std::pair<double, double> origin{
/*
std::pair<double, double> utmToLatitudeLongitude(
     const std::pair<double, double> &xy,
     int utmProjectionZone = 10);
*/

}
}
