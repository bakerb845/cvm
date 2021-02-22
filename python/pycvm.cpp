#include "cvm/version.hpp"
#include "cvm/enums.hpp"
#include "include/pConstants.hpp"
#include "include/pGeodetic.hpp"
#include <pybind11/pybind11.h>

PYBIND11_MODULE(pycvm, m)
{
    m.attr("__version__") = CVM_VERSION;   
    // Layer identifiers
    pybind11::enum_<CVM::LayerIdentifier> (m, "LayerIdentifier")
        .value("top", CVM::LayerIdentifier::TOP,
               "The upper crustal layer of the CVM")
        .value("middle", CVM::LayerIdentifier::MIDDLE,
               "The middle-to-lower crustal layer of the CVM")
        .value("bottom", CVM::LayerIdentifier::BOTTOM,
               "The mantle layer of the CVM");
    // Output model type
    pybind11::enum_<CVM::FileType> (m, "FileType")
        .value("nll", CVM::FileType::NLL,
               "NonLinLoc binary file format")
        .value("vtk", CVM::FileType::VTK,
               "Legacy VTK format");
    // Model constants
    PCVM::initializeConstants(m);
    // Geodetic calculations
    PCVM::initializeGeodetic(m);
}
