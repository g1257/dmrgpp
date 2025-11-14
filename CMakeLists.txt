cmake_minimum_required(VERSION 3.28)

project(
  dmrgpp
  VERSION 1.1
  LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${dmrgpp_SOURCE_DIR}/CMake)

find_package(MPI COMPONENTS CXX)
include(CTest)

find_package(MPI COMPONENTS CXX)
include(HandleHDF5)
find_package(BLAS)
find_package(LAPACK)

add_subdirectory(external/Catch2)
include(Catch)

add_subdirectory(PsimagLite)
add_subdirectory(src)
# add_subdirectory(doc)
