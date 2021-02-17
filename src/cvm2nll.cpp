#include <iostream>
#include <string>
#include <getopt.h>
#include "cvm/model.hpp"
#include "cvm/options.hpp"
#include "cvm/selection.hpp"

/*
struct ProgramOptions
{
    double latMin = 46.75;
    double latMax =-123.4;
    double lonMin = 48.25;
    double lonMax =-121.4;
    double depMin = 0;
    double depMax = 20;
};
*/

int main(int argc, char *argv[])
{
    // Parse command line arguments
    const struct option longOptions[] = 
    {
        {"ini",  required_argument, 0, 'i'},
        {"help", no_argument,       0, '?'},
        {"help", no_argument,       0, 'h'},
        {0,      0,                 0, 0},
    };
    std::string optionFile = "";
    int optionIndex = 0;
    while (true)
    {
        auto c = getopt_long(argc, argv, "i:?h", longOptions, &optionIndex);
        if (c == 'h' || c == '?')
        {
            std::cout << "Usage:" << std::endl;
            std::cout << "cvm2nll --ini fileName" << std::endl;
            std::cout << "  --ini name of initialization file" << std::endl;
            std::cout << "  --help prints this message" << std::endl;
            return EXIT_SUCCESS;
        }
        else if (c == 'i')
        {
            optionFile = std::string(optarg); //longOptions[optionIndex].name);
            break; 
        }
        else
        {
            std::cerr << "Usage:" << std::endl;
            std::cerr << "cvm2nll --ini fileName" << std::endl;
            std::cerr << "  --ini name of initialization file" << std::endl;
            std::cerr << "  --help prints this message" << std::endl;
            return EXIT_FAILURE;
        }
    }
    // Parse the ini file
    std::cout << "Parsing: " << optionFile << std::endl;
    CVM::Options options;
    try
    {
        options.parse(optionFile);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
 
    // Load the file
    CVM::Model model;
    model.load(options);
/*
    auto selection = options.getSelection();
std::string fileNameLayer1 = "../data/vp_16Bl1.bin";
std::string fileNameLayer2 = "../data/vp_16Bl2.bin";
std::string fileNameLayer3 = "../data/vp_16Bl3.bin";
    std::cout << "Loading CVM..." << std::endl;
    CVM::Layer<CVM::LayerIdentifier::TOP> pLayer1;
    CVM::Layer<CVM::LayerIdentifier::MIDDLE> pLayer2;
    CVM::Layer<CVM::LayerIdentifier::BOTTOM> pLayer3;
    pLayer1.load(fileNameLayer1, selection, true);
    pLayer2.load(fileNameLayer2, selection, true);
    pLayer3.load(fileNameLayer3, selection, true);
    std::cout << "Merging layers and interpolating..." << std::endl; 

    std::cout << "Writing to NLL format..." << std::endl;
*/
    return EXIT_SUCCESS;
}

