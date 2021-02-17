#include <iostream>
#include <vector>
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
                      const Layer<CVM::LayerIdentifier::BOTTOM> &layer3)
{
    const float *vt = nullptr;
    vt = layer1.getPVelocityPointer(); 
    if (lvp){vt = layer1.getSVelocityPointer();}

    const float *vm = nullptr;
    vm = layer2.getPVelocityPointer(); 
    if (lvp){vm = layer2.getSVelocityPointer();}

    const float *vb = nullptr;
    vb = layer3.getPVelocityPointer(); 
    if (lvp){vb = layer3.getSVelocityPointer();}

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
    layer1.load(options.getPVelocityFileName(LayerIdentifier::TOP),
                selection, isP);
    layer2.load(options.getPVelocityFileName(LayerIdentifier::MIDDLE),
                selection, isP);
    layer3.load(options.getPVelocityFileName(LayerIdentifier::BOTTOM),
                selection, isP);
    // Heal the layer
    std::cout << "Healing P velocity model..." << std::endl;

    // Release memory
    layer1.clear();
    layer2.clear();
    layer3.clear();
    // S velocity model
    isP = false; 
    layer1.load(options.getSVelocityFileName(LayerIdentifier::TOP),
                selection, isP);
    layer2.load(options.getSVelocityFileName(LayerIdentifier::MIDDLE),
                selection, isP);
    layer3.load(options.getSVelocityFileName(LayerIdentifier::BOTTOM),
                selection, isP);


}
