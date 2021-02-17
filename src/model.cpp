#include <iostream>
#include <vector>
#include "cvm/model.hpp"
#include "cvm/selection.hpp"
#include "cvm/options.hpp"
#include "cvm/layer.hpp"

using namespace CVM;

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
