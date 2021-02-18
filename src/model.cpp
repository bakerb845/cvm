#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include <cassert>
#include "cvm/model.hpp"
#include "cvm/selection.hpp"
#include "cvm/constants.hpp"
#include "cvm/options.hpp"
#include "cvm/layer.hpp"
#include "private/interpolate.hpp"

using namespace CVM;

namespace
{

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

    auto dx2 = Constants::getGridSpacingInX(LayerIdentifier::MIDDLE);
    auto dy2 = Constants::getGridSpacingInY(LayerIdentifier::MIDDLE);
    auto dz2 = Constants::getGridSpacingInZ(LayerIdentifier::MIDDLE);
    auto nx2 = layer2.getNumberOfGridPointsInX(); 
    auto ny2 = layer2.getNumberOfGridPointsInY();
    auto nz2 = layer2.getNumberOfGridPointsInZ();

    auto dx3 = Constants::getGridSpacingInX(LayerIdentifier::BOTTOM);
    auto dy3 = Constants::getGridSpacingInY(LayerIdentifier::BOTTOM);
    auto dz3 = Constants::getGridSpacingInZ(LayerIdentifier::BOTTOM);
    auto nx3 = layer3.getNumberOfGridPointsInX(); 
    auto ny3 = layer3.getNumberOfGridPointsInY();
    auto nz3 = layer3.getNumberOfGridPointsInZ();

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
                auto x = static_cast<float> (dx*ix);
                auto y = static_cast<float> (dy*iy);
                auto z = static_cast<float> (z0 + dz*iz);
                // Check if z is in one of the problem areas

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
    std::cout << "Healing P velocity model..." << std::endl;
    double z0 = 0;
    double z1 = 40*1000;
    double dx = options.getNLLGridSpacingInX();
    double dy = options.getNLLGridSpacingInY();
    double dz = options.getNLLGridSpacingInZ();
    int nx, ny, nz;
    interpolateModel(isP,
                     z0, z1, 
                     dx, dy, dz,
                     layer1, layer2, layer3,
                     &nx, &ny, &nz, &pImpl->mPModel);
    // Release memory
    layer1.clear();
    layer2.clear();
    layer3.clear();
    // S velocity model
    isP = false; 
    layer1.load(options, isP);
    layer2.load(options, isP);
    layer3.load(options, isP);


}
