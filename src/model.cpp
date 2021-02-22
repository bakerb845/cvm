#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <limits>
#include <cassert>
#include "cvm/model.hpp"
#include "cvm/selection.hpp"
#include "cvm/constants.hpp"
#include "cvm/options.hpp"
#include "cvm/layer.hpp"
#include "private/interpolate.hpp"
#include "../external/visit_writer.h"

using namespace CVM;

namespace
{

float interp3d(const int nxLayer, const int nyLayer, const int nzLayer,
               const double dxLayer,
               const double dyLayer,
               const double dzLayer,
               const float x, const float y, const float z,
               const float *vLayer)
{
    auto ixq = static_cast<int> (x/dxLayer);
    auto iyq = static_cast<int> (y/dyLayer);
    auto izq = static_cast<int> (z/dzLayer);
    if (ixq == nxLayer - 1){ixq = ixq - 1;}
    if (iyq == nyLayer - 1){iyq = iyq - 1;} 
    if (izq == nzLayer - 1){izq = izq - 1;} 
#ifndef NDEBUG
    assert(ixq >= 0 && ixq < nxLayer - 1);
    assert(iyq >= 0 && iyq < nyLayer - 1);
    assert(izq >= 0 && izq < nzLayer - 1);
#endif
    auto indx = izq*nxLayer*nyLayer + iyq*nxLayer + ixq;
    std::array<float, 8> v8;
    v8[0] = vLayer[indx];
    v8[1] = vLayer[indx + 1]; 
    v8[2] = vLayer[indx + nxLayer];
    v8[3] = vLayer[indx + nxLayer + 1]; 

    v8[4] = vLayer[indx + nxLayer*nyLayer];
    v8[5] = vLayer[indx + nxLayer*nyLayer + 1]; 
    v8[6] = vLayer[indx + nxLayer*nyLayer + nxLayer];
    v8[7] = vLayer[indx + nxLayer*nyLayer + nxLayer + 1]; 
 
    auto vInt = trilinearInterpolation(static_cast<float> (ixq*dxLayer),
                                       static_cast<float> ((ixq + 1)*dxLayer),
                                       static_cast<float> (iyq*dyLayer),
                                       static_cast<float> ((iyq + 1)*dyLayer),
                                       static_cast<float> (izq*dzLayer),
                                       static_cast<float> ((izq + 1)*dzLayer),
                                       x, y, z,
                                       v8.data());
    return vInt;

}

void interpolateModel(const bool lvp,
                      const double z0, const double z1,
                      const double dx, const double dy, const double dz,
                      const Layer<CVM::LayerIdentifier::TOP> &layer1,
                      const Layer<CVM::LayerIdentifier::MIDDLE> &layer2,
                      const Layer<CVM::LayerIdentifier::BOTTOM> &layer3,
                      int *nx, int *ny, int *nz,
                      std::vector<float> *v)
{
    const float *vt = nullptr;
    const float *vm = nullptr;
    const float *vb = nullptr;
    if (lvp)
    {
        vt = layer1.getPVelocityPointer(); 
        vm = layer2.getPVelocityPointer(); 
        vb = layer3.getPVelocityPointer();
    }
    else
    {
        vt = layer1.getSVelocityPointer();
        vm = layer2.getSVelocityPointer();
        vb = layer3.getSVelocityPointer();
    }

    auto [x0, y0] = layer1.getSouthWestCornerUTM();
    auto [x1, y1] = layer1.getNorthEastCornerUTM();
#ifndef NDEBUG
    auto [x0Test2, y0Test2] = layer2.getSouthWestCornerUTM();
    auto [x1Test2, y1Test2] = layer2.getNorthEastCornerUTM();
    assert(std::abs(x0 - x0Test2) < 1.e-8);
    assert(std::abs(y0 - y0Test2) < 1.e-8);
    assert(std::abs(x1 - x1Test2) < 1.e-8);
    assert(std::abs(y1 - y1Test2) < 1.e-8);
    auto [x0Test3, y0Test3] = layer3.getSouthWestCornerUTM();
    auto [x1Test3, y1Test3] = layer3.getNorthEastCornerUTM();
    assert(std::abs(x0 - x0Test3) < 1.e-8);
    assert(std::abs(y0 - y0Test3) < 1.e-8);
    assert(std::abs(x1 - x1Test3) < 1.e-8);
    assert(std::abs(y1 - y1Test3) < 1.e-8);
#endif
    // Figure out grid sizes
    auto dx1 = Constants::getGridSpacingInX(LayerIdentifier::TOP);
    auto dy1 = Constants::getGridSpacingInY(LayerIdentifier::TOP);
    auto dz1 = Constants::getGridSpacingInZ(LayerIdentifier::TOP);
    auto nx1 = layer1.getNumberOfGridPointsInX(); 
    auto ny1 = layer1.getNumberOfGridPointsInY();
    auto nz1 = layer1.getNumberOfGridPointsInZ();
    auto zStart1 = Constants::getLayerStartDepth(LayerIdentifier::TOP);
    auto zEnd1   = Constants::getLayerEndDepth(LayerIdentifier::TOP); 

    auto dx2 = Constants::getGridSpacingInX(LayerIdentifier::MIDDLE);
    auto dy2 = Constants::getGridSpacingInY(LayerIdentifier::MIDDLE);
    auto dz2 = Constants::getGridSpacingInZ(LayerIdentifier::MIDDLE);
    auto nx2 = layer2.getNumberOfGridPointsInX(); 
    auto ny2 = layer2.getNumberOfGridPointsInY();
    auto nz2 = layer2.getNumberOfGridPointsInZ();
    auto zStart2 = Constants::getLayerStartDepth(LayerIdentifier::MIDDLE);
    auto zEnd2   = Constants::getLayerEndDepth(LayerIdentifier::MIDDLE); 


    auto dx3 = Constants::getGridSpacingInX(LayerIdentifier::BOTTOM);
    auto dy3 = Constants::getGridSpacingInY(LayerIdentifier::BOTTOM);
    auto dz3 = Constants::getGridSpacingInZ(LayerIdentifier::BOTTOM);
    auto nx3 = layer3.getNumberOfGridPointsInX(); 
    auto ny3 = layer3.getNumberOfGridPointsInY();
    auto nz3 = layer3.getNumberOfGridPointsInZ();
    auto zStart3 = Constants::getLayerStartDepth(LayerIdentifier::BOTTOM);
    auto zEnd3   = Constants::getLayerEndDepth(LayerIdentifier::BOTTOM); 

    // Try to pack as many x and y grid points in the available space
    *nx = 0;
    for (int i = 1; i < std::numeric_limits<int>::max(); ++i)
    {
        if (i*dx > (nx1 - 1)*dx1)
        {
            *nx = i;
            break;
        }
    }
    *ny = 0;
    for (int i = 1; i < std::numeric_limits<int>::max(); ++i)
    {
        if (i*dy > (ny1 - 1)*dy1)
        {
            *ny = i;
            break;
        }
    }
#ifndef NDEBUG
    assert(*nx > 0);
    assert(*ny > 0);
#endif
    *nz = static_cast<int> (std::round(z1 - z0)/dz) + 1;
    v->resize((*nx)*(*ny)*(*nz), 0);
    // Begin interpolation process
    for (int iz = 0; iz < *nz; ++iz)
    {
        for (int iy = 0; iy < *ny; ++iy)
        {
            for (int ix = 0; ix < *nx; ++ix)
            {
                auto idst = iz*(*nx)*(*ny) + iy*(*nx) + ix;
                auto x = static_cast<float> (dx*ix);
                auto y = static_cast<float> (dy*iy);
                auto z = static_cast<float> (z0 + dz*iz);
                // Check if z is in a layer
                if (z >= zStart1 && z <= zEnd1)
                {
                    v->at(idst) = interp3d(nx1, ny1, nz1,
                                           dx1, dy1, dz1,
                                           x, y,
                                           z - static_cast<float> (zStart1),
                                           vt);
                }
                else if (z >= zStart2 && z <= zEnd2)
                {
                    v->at(idst) = interp3d(nx2, ny2, nz2,
                                           dx2, dy2, dz2,
                                           x, y,
                                           z - static_cast<float> (zStart2),
                                           vm);
                }
                else if (z >= zStart3 && z <= zEnd3)
                {
                    v->at(idst) = interp3d(nx3, ny3, nz3,
                                           dx3, dy3, dz3,
                                           x, y,
                                           z - static_cast<float> (zStart3),
                                           vt);
                }
                else
                {

                }
            }
        }
    }
}

}

class Model::ModelImpl
{
public:
    std::vector<float> mPModel;
    std::vector<float> mSModel;
    double mDx = 0;
    double mDy = 0;
    double mDz = 0;
    double mX0 = 0;
    double mY0 = 0;
    double mZ0 = 0;
    int mX = 0;
    int mY = 0;
    int mZ = 0;
    bool mLoaded = false;
};

/// C'tor
Model::Model() :
    pImpl(std::make_unique<ModelImpl>())
{
}

/// Copy c'tor
Model::Model(const Model &model)
{
    *this = model;
}

/// Move c'tor
Model::Model(Model &&model) noexcept
{
    *this = std::move(model);
}

/// Copy assignment
Model& Model::operator=(const Model &model)
{
    if (&model == this){return *this;}
    pImpl = std::make_unique<ModelImpl> (*model.pImpl);
    return *this;
}

/// Move assignment
Model& Model::operator=(Model &&model) noexcept
{
    if (&model == this){return *this;}
    pImpl = std::move(model.pImpl);
    return *this;
}

/// Clear model and release memory
void Model::clear() noexcept
{
    pImpl->mPModel.clear();
    pImpl->mSModel.clear();
    pImpl->mDx = 0;
    pImpl->mDy = 0;
    pImpl->mDz = 0;
    pImpl->mX0 = 0;
    pImpl->mY0 = 0;
    pImpl->mZ0 = 0;
    pImpl->mX = 0;
    pImpl->mY = 0;
    pImpl->mZ = 0;
    pImpl->mLoaded = false;
}

/// Destructor
Model::~Model() = default;

/// Loads the CVM
void Model::load(const Options &options)
{
    auto selection = options.getSelection();
    // Load each layer of the model 
    std::cout << "Loading P velocity model..." << std::endl;
    Layer<CVM::LayerIdentifier::TOP> layer1;
    Layer<CVM::LayerIdentifier::MIDDLE> layer2;
    Layer<CVM::LayerIdentifier::BOTTOM> layer3;
    // P velocity model first
    bool isP = true;
    layer1.load(options, isP);
    layer2.load(options, isP);
    layer3.load(options, isP);
    // Heal the layer
    std::cout << "Interpolating P velocity model onto regular grid..."
              << std::endl;
    auto zStartStop = selection.getMinimumAndMaximumDepth();
    double z0 = zStartStop.first;
    double z1 = zStartStop.second;
std::cout << z0 << " " << z1 << std::endl;
    double dx = options.getNLLGridSpacingInX();
    double dy = options.getNLLGridSpacingInY();
    double dz = options.getNLLGridSpacingInZ();
    int nx, ny, nz;
    interpolateModel(isP,
                     z0, z1, 
                     dx, dy, dz,
                     layer1, layer2, layer3,
                     &nx, &ny, &nz, &pImpl->mPModel);
    // S velocity model
    std::cout << "Loading S velocity model..." << std::endl;
    isP = false; 
    layer1.load(options, isP);
    layer2.load(options, isP);
    layer3.load(options, isP);
    std::cout << "Interpolating S velocity model onto regular grid..."
              << std::endl;
    interpolateModel(isP,
                     z0, z1,
                     dx, dy, dz,
                     layer1, layer2, layer3,
                     &nx, &ny, &nz, &pImpl->mSModel);

    // Update model information
    auto origin = layer1.getSouthWestCornerUTM();
    pImpl->mX0 = origin.first;
    pImpl->mY0 = origin.second;
    pImpl->mZ0 = z0;
    pImpl->mDx = dx;
    pImpl->mDy = dy;
    pImpl->mDz = dz;
    pImpl->mX = nx;
    pImpl->mY = ny;
    pImpl->mZ = nz;
    // Release memory
    layer1.clear();
    layer2.clear();
    layer3.clear();

    pImpl->mLoaded = true;
}

/// Is the model loaded
bool Model::isLoaded() const noexcept
{
    return pImpl->mLoaded;
}

/// Number of grid points
int Model::getNumberOfGridPointsInX() const
{
    if (!isLoaded()){throw std::runtime_error("Model not loaded");}
    return pImpl->mX;
}

int Model::getNumberOfGridPointsInY() const
{
    if (!isLoaded()){throw std::runtime_error("Model not loaded");}
    return pImpl->mY;
}

int Model::getNumberOfGridPointsInZ() const
{
    if (!isLoaded()){throw std::runtime_error("Model not loaded");}
    return pImpl->mZ;
}

/// Grid spacing
double Model::getGridSpacingInX() const
{
    if (!isLoaded()){throw std::runtime_error("Model not loaded");}
    return pImpl->mDx;
}

double Model::getGridSpacingInY() const
{
    if (!isLoaded()){throw std::runtime_error("Model not loaded");}
    return pImpl->mDy;
}

double Model::getGridSpacingInZ() const
{
    if (!isLoaded()){throw std::runtime_error("Model not loaded");}
    return pImpl->mDz;
}

const float *Model::getPVelocityPointer() const
{
    if (!isLoaded()){throw std::runtime_error("Model not loaded");}
    return pImpl->mPModel.data();
}

const float *Model::getSVelocityPointer() const
{
    if (!isLoaded()){throw std::runtime_error("Model not loaded");}
    return pImpl->mSModel.data();
}

void Model::writePVelocitiesToFile(const std::string &fileName,
                                   const FileType fileType)
{
    auto nx = getNumberOfGridPointsInX();
    auto ny = getNumberOfGridPointsInY();
    auto nz = getNumberOfGridPointsInZ();
    auto dx = getGridSpacingInX();
    auto dy = getGridSpacingInY();
    auto dz = getGridSpacingInZ();
    auto vp = getPVelocityPointer(); 
    auto vs = getSVelocityPointer();
    if (fileType == FileType::NLL)
    {
    }
    else if (fileType == FileType::VTK)
    {
        const int useBinary = static_cast<int> (true);
        const int nVars = 2;
        std::array<int, 3> dims{nz, ny, nx};
        const char *const varNames[2] = {"vp", "vs"};
        const float *vars[2] = {vp, vs}; 
        std::array<int, 2> centering{0, 0};
        std::array<int, 2> varDim{1, 1};
        std::vector<float> x(nx, 0);
        std::vector<float> y(ny, 0);
        std::vector<float> z(nz, 0);
        for (int i = 0; i < nx; ++i){x[i] = i*dx;}
        for (int i = 0; i < ny; ++i){y[i] = i*dy;}
        for (int i = 0; i < nz; ++i){z[i] = i*dz;} 
        write_rectilinear_mesh("test", useBinary,
                               dims.data(), z.data(), y.data(), x.data(),
                               nVars, varDim.data(), centering.data(),
                               varNames, vars);
    }
}
