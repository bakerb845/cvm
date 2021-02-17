#include <iostream>
#include <cmath>
#include <vector>
#include "private/interpolate.hpp"
#include <gtest/gtest.h>
namespace
{

template<typename T>
T fxy(const T x, const T y)
{
    return 2*x + static_cast<T> (0.1)*y;
}

template<typename T>
T fxyz(const T x, const T y, const T z)
{
    return 3*x + static_cast<T> (0.2)*y + static_cast<T> (-0.4)*z;
}

template<typename T>
void makeGrid(const int nx, const int ny, const int nz,
              const T x0, const T dx, 
              const T y0, const T dy, 
              const T z0, const T dz,
              std::vector<T> &x, 
              std::vector<T> &y,
              std::vector<T> &z,
              std::vector<T> &f,
              const bool phase1 = true)
{
    int nxyz = nx*ny*nz;
    if (phase1)
    {
        x.resize(nx, 0);
        y.resize(ny, 0);
        z.resize(nz, 0);
    }
    else
    {
        x.resize(nxyz, 0); 
        y.resize(nxyz, 0); 
        z.resize(nxyz, 0); 
    }
    f.resize(nxyz, 0);
    for (int iz=0; iz<nz; ++iz)
    {
        for (int iy=0; iy<ny; ++iy)
        {
            for (int ix=0; ix<nx; ++ix)
            {
                auto indx = iz*nx*ny + iy*nx + ix;
                auto ixuse = indx;
                auto iyuse = indx;
                auto izuse = indx;
                if (phase1 == true)
                {
                    ixuse = ix;
                    iyuse = iy;
                    izuse = iz;
                }
                x[ixuse] = x0 + ix*dx;
                y[iyuse] = y0 + iy*dy;
                z[izuse] = z0 + iz*dz;
                f[indx] = fxyz(x[ixuse], y[iyuse], z[izuse]);
            }
        }
    }
}

template<typename T>
void makeGrid(const int nx, const int ny, 
              const T x0, const T dx, 
              const T y0, const T dy, 
              std::vector<T> &x, 
              std::vector<T> &y, 
              std::vector<T> &f, 
              bool phase1 = true)
{
    auto nxy = nx*ny;
    if (phase1)
    {
        x.resize(nx, 0); 
        y.resize(ny, 0); 
    }   
    else
    {
        x.resize(nxy, 0); 
        y.resize(nxy, 0); 
    }   
    f.resize(nxy, 0); 
    for (int iy=0; iy<ny; ++iy)
    {
        for (int ix=0; ix<nx; ++ix)
        {
            auto indx = iy*nx + ix; 
            auto ixuse = indx;
            auto iyuse = indx;
            if (phase1)
            {
                ixuse = ix; 
                iyuse = iy; 
            }
            x[ixuse] = x0 + ix*dx;
            y[iyuse] = y0 + iy*dy;
            f[indx] = fxy(x[ixuse], y[iyuse]);
        }   
    }   
}

TEST(Interpolate, Interpolate2D)
{
    int nx = 5;
    int ny = 6;
    double dx = 0.2;
    double dy = 0.5;
    double x0 = 1;
    double y0 = 2;
    double x1 = x0 + (nx - 1)*dx;
    double y1 = y0 + (ny - 1)*dy;
    int nxInt = 15; 
    int nyInt = 18;
    double dxInt = (x1 - x0)/(nxInt - 1);
    double dyInt = (y1 - y0)/(nyInt - 1);
    std::vector<double> x, y, f;
    makeGrid(nx, ny, x0, dx, y0, dy, x, y, f, true);
    std::vector<double> xq, yq, fqRef;
    makeGrid(nxInt, nyInt, x0, dxInt, y0, dyInt, xq, yq, fqRef, true);
    double v[4];
    double resMax = 0;
    for (int iy = 0; iy < nyInt; ++iy)
    {
        for (int ix = 0; ix < nxInt; ++ix)
        {
            auto ixq = static_cast<int> ((xq[ix] - x0)/dx);
            auto iyq = static_cast<int> ((yq[iy] - y0)/dy);
            if (ixq == nx - 1){ixq = ixq - 1;}
            if (iyq == ny - 1){iyq = iyq - 1;}
            auto indx = iyq*nx + ixq;
            v[0] = f.at(indx);
            v[1] = f.at(indx + 1);
            v[2] = f.at(indx + nx);
            v[3] = f.at(indx + nx + 1);
            auto xInt1 = x[ixq]; //x0 + ixq*dx;
            auto yInt1 = y[iyq]; //y0 + iyq*dy;
            auto xInt2 = x[ixq + 1]; //x0 + (ixq + 1)*dx;
            auto yInt2 = y[iyq + 1]; //y0 + (iyq + 1)*dy;
            auto fq = bilinearInterpolation(xInt1, xInt2,
                                            yInt1, yInt2,
                                            xq[ix], yq[iy], v);
            auto res = fqRef.at(iy*nxInt + ix) - fq;
            resMax = std::max(std::abs(res), resMax);
        }
    }
    EXPECT_NEAR(resMax, 0, 1.e-14);
    //std::cout << resMax << std::endl; 
}

TEST(Interpolate, Interpolate3D)
{
    int nx = 5;
    int ny = 6;
    int nz = 7;
    double dx = 0.2;
    double dy = 0.5;
    double dz = 0.7;
    double x0 = 0;
    double y0 = 0;
    double z0 = 0;
    double x1 = x0 + (nx - 1)*dx;
    double y1 = y0 + (ny - 1)*dy;
    double z1 = z0 + (nz - 1)*dz;
    int nxInt = 15; 
    int nyInt = 18; 
    int nzInt = 21;
    double dxInt = (x1 - x0)/(nxInt - 1); 
    double dyInt = (y1 - y0)/(nyInt - 1); 
    double dzInt = (z1 - z0)/(nzInt - 1);
    std::vector<double> x, y, z, f;
    makeGrid(nx, ny, nz,
             x0, dx, y0, dy, z0, dz,
             x, y, z, f, true);
    std::vector<double> xq, yq, zq, fqRef;
    makeGrid(nxInt, nyInt, nzInt,
             x0, dxInt, y0, dyInt, z0, dzInt,
             xq, yq, zq, fqRef, true);
    double v[8];
    double resMax = 0;
    for (int iz = 0; iz < nzInt; ++iz)
    {
        for (int iy = 0; iy < nyInt; ++iy)
        {   
            for (int ix = 0; ix < nxInt; ++ix)
            {
                auto ixq = static_cast<int> ((xq[ix] - x0)/dx);
                auto iyq = static_cast<int> ((yq[iy] - y0)/dy);
                auto izq = static_cast<int> ((zq[iz] - z0)/dz);
                if (ixq == nx - 1){ixq = ixq - 1;} 
                if (iyq == ny - 1){iyq = iyq - 1;} 
                if (izq == nz - 1){izq = izq - 1;}
                auto indx = izq*nx*ny + iyq*nx + ixq;
                v[0] = f.at(indx);
                v[1] = f.at(indx + 1); 
                v[2] = f.at(indx + nx);
                v[3] = f.at(indx + nx + 1); 

                v[4] = f.at(indx + nx*ny);
                v[5] = f.at(indx + nx*ny + 1);
                v[6] = f.at(indx + nx*ny + nx);
                v[7] = f.at(indx + nx*ny + nx + 1);
                auto xInt0 = x[ixq];
                auto yInt0 = y[iyq];
                auto zInt0 = z[izq];
                auto xInt1 = x[ixq + 1];
                auto yInt1 = y[iyq + 1];
                auto zInt1 = z[izq + 1];
                auto fq = trilinearInterpolation(xInt0, xInt1,
                                                 yInt0, yInt1,
                                                 zInt0, zInt1,
                                                 xq[ix], yq[iy], zq[iz],
                                                 v);
                auto res = fqRef.at(iz*nxInt*nyInt + iy*nxInt + ix) - fq;
                resMax = std::max(std::abs(res), resMax);
            }
        }
    }
    EXPECT_NEAR(resMax, 0, 1.e-14);
    //std::cout << resMax << std::endl;
}
}
