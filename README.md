# Overview 

This is a utility for converting the Cascadia Community Velocity Model to a NonLinLoc file.

# Prerequisites

    1.  A C++17 compiler.
    2.  [CMake](https://cmake.org/) v3.10 or greater.
    3.  [Boost](https://www.boost.org/)
    3.  A lot of RAM.  Worst case is this program will about 10-15 Gb of RAM.  It's a big model.
    4.  [pybind11](https://github.com/pybind/pybind11) and Python3 for creating a Python interface.

# Configure 

An example configuration script is given as config.sh.  Effectively all it does is create a build directory and instantiate CMake.

An example configuration script, config.sh, would look like

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


