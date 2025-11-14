cmake_minimum_required(VERSION 3.31)

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

include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.11.0
  FIND_PACKAGE_ARGS NAMES Catch2)
FetchContent_MakeAvailable(Catch2)
include(Catch)

add_subdirectory(PsimagLite)
add_subdirectory(src)
# add_subdirectory(doc)
