# Overview 

This is a utility for converting the Cascadia Community Velocity Model to a NonLinLoc file.

# Prerequisites

    1.  A C++17 compiler.
    2.  [CMake](https://cmake.org/) v3.10 or greater.
    3.  [Boost](https://www.boost.org/)
    3.  A lot of RAM.  Worst case is this program will about 10-15 Gb of RAM.  It's a big model.
    4.  Optionally, [pybind11](https://github.com/pybind/pybind11) and Python3 for creating a Python interface.

# Download the Software

Descend into the working directory of your choice, clone the repository, then descend into the CVM directory

    git clone https://github.com/bakerb845/cvm.git
    cd cvm

# Configure 

In theory, CMake can automatically create a Makefile for your machine.  In practice, this is rarely realized, especially if  things are not installed in standard locations.  If you'd like to roll the dice then try

    mkdir build
    cd build
    cmake ..

An example configuration script, config.sh, for building the entire may look like

    #!/usr/bin/bash
    BUILD_DIR=clang_build
    if [ -d ${BUILD_DIR} ]; then
       rm -rf ${BUILD_DIR}
    fi
    mkdir ${BUILD_DIR}
    cd ${BUILD_DIR}
    cmake .. \
    -DCMAKE_CXX_COMPILER=clang++ \
    -DCMAKE_C_COMPILER=clang \
    -DCMAKE_CXX_FLAGS="-Wall -O2" \
    -DCMAKE_C_COMPILER="-Wall -O2" \
    -DPYTHON_EXECUTABLE="/home/bbaker/anaconda3/bin/python" \
    -DPYTHON_LIBRARY="/home/bbaker/anaconda3/lib/libpython3.8.so"

If you do not want the Python bindings then you can disregard the given Python information.

# Building the Software

After configuring the software, descend into the build directory and compile.

    make
    make test
    make install

The last command may require sudo privileges.

# Usage

## Command Line Interface

One way to run the software is through the command line interface.  This requires editing an initialization file.  An example ini file is provided in examples/seattle.ini.  To use this you could copy this example ini file to your local directory, modify the file paths and options as required, and extract the CVM with something analogous to

    cp ROOT_SOURCE_DIRECTORY/examples/seattle.ini ./seattle_modify.ini
    cvm2nll --ini seattle_modify.ini 

## Python

If you have compiled the Python bindings and set the PYTHONPATH in your bashrc then you can leverage some library functionality from Python.  This is useful for when you, say, have station latitudes and longitudes that you must then make consistent with a Cartesian system in NonLinLoc.  unit\_test.py will provide examples of how to import the library and use different modules.

    
