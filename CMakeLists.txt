cmake_minimum_required(VERSION 3.28)

project(
  dmrgpp
  VERSION 1.1
  LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${dmrgpp_SOURCE_DIR}/CMake)

set(QMC_CXX_STANDARD
    11
    CACHE STRING "DMRG++ C++ language standard")
if(NOT QMC_CXX_STANDARD EQUAL 11)
  message(
    WARNING
      "C++11 is the only language standard officially supported by dmrpp. "
      "Using other versions of the C++ standard is unsupported and done entirely at user's own risk."
  )
endif()

include(HandleHDF5)

find_package(BLAS)
find_package(LAPACK)

add_subdirectory(src)
# add_subdirectory(doc)
