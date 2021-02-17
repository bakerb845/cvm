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

    ./config.sh


