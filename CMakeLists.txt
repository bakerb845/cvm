cmake_minimum_required(VERSION 3.10)
project(cvm VERSION 0.1.0 LANGUAGES CXX C)
enable_testing()


include(CheckCXXCompilerFlag)
find_package(Boost)
find_package(GTest)

# Do I have std::filesystem?
try_compile(HAVE_STD_FILESYSTEM ${CMAKE_BINARY_DIR}/temp ${CMAKE_SOURCE_DIR}/testing/fs.cpp
            LINK_LIBRARIES -lstdc++fs
            CXX_STANDARD 17)
if (${HAVE_STD_FILESYSTEM})
   message("Using std::filesystem")
   add_compile_definitions(USE_STD_FILESYSTEM)
   link_libraries(stdc++fs)
else()
   message("Disabling filesystem") 
endif()

configure_file(${CMAKE_SOURCE_DIR}/include/private/version.hpp.in
               ${CMAKE_SOURCE_DIR}/include/cvm/version.hpp)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${Boost_INCLUDE_DIRS}
)

set(CVM_LIBRARY_SRC
    src/constants.cpp
    src/layer.cpp
    src/geodetic.cpp
    src/selection.cpp
    src/model.cpp
    src/version.cpp
    external/visit_writer.c)
if (Boost_FOUND)
   set(CVM_LIBRARY_SRC ${CVM_LIBRARY_SRC} src/options.cpp)
endif()
add_library(cvm SHARED ${CVM_LIBRARY_SRC})
set_target_properties(cvm PROPERTIES
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)

set(CVM2NLL_SRC
    src/cvm2nll.cpp)
add_executable(cvm2nll ${CVM2NLL_SRC})
set_target_properties(cvm2nll PROPERTIES
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)
target_link_libraries(cvm2nll PRIVATE cvm)


if (${GTEST_FOUND})
   message("Gtest found - will compile unit tests")
   set(TEST_SRC
       testing/main.cpp
       testing/geodetic.cpp
       testing/constants.cpp
       testing/interpolate.cpp)
   add_executable(unitTests ${TEST_SRC})
   set_target_properties(unitTests PROPERTIES
                         CXX_STANDARD 17
                         CXX_STANDARD_REQUIRED 17
                         CXX_EXTENSIONS NO)
   target_link_libraries(unitTests cvm ${GTEST_BOTH_LIBRARIES})
   add_test(NAME unitTests COMMAND unitTests)
endif()

if (${BUILD_PYTHON})
   find_package(pybind11 REQUIRED)
   message("Will create Python interface")
   file(COPY ${CMAKE_SOURCE_DIR}/python/unit_test.py DESTINATION .)
   #find_package(PythonInterp 3)
   #find_package(PythonLibs 3)
   add_library(pycvm MODULE
               python/pycvm.cpp
               python/constants.cpp
               python/geodetic.cpp)
   target_link_libraries(pycvm PRIVATE pybind11::module cvm) 
   set_target_properties(pycvm PROPERTIES
                         CXX_STANDARD 17
                         CXX_STANDARD_REQUIRED 17
                         CXX_EXTENSIONS NO
                         PREFIX "")
   target_include_directories(pycvm PRIVATE pybind11
                              ${CMAKE_SOURCE_DIR}/python
                              ${CMAKE_SOURCE_DIR}/python/include)
endif()

################################################################################
#                                Installation                                  #
################################################################################
include(GNUInstallDirs)
if (${BUILD_PYTHON})
   install(TARGETS cvm2nll cvm pycvm
           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
           LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
           ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
           PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
else()
   install(TARGETS cvm2nll cvm
           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
           LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
           ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
           PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

##########################################################################################
#                                     CPACK Packaging                                    #
##########################################################################################
set(CPACK_PACKAGE_NAME "cvm")
set(CPACK_PACKAGE_VENDOR "UUSS")
set(CPACK_PACKAGE_CONTACT "bbaker@seis.utah.edu")
set(CPACK_PACKAGE_LICENSE "MIT")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Tools for unpacking the Cascadia Community Velocity Model")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
