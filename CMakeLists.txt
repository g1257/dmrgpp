cmake_minimum_required(VERSION 3.31)

project(
  dmrgpp
  VERSION 1.1
  LANGUAGES CXX)

# Set default to RelWithDebInfo if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      "RelWithDebInfo"
      CACHE STRING "Choose the type of build" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

option(BUILD_TESTING "Enable building tests" ON)
if(BUILD_TESTING)
  enable_testing()

  include(FetchContent)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.11.0
    SYSTEM FIND_PACKAGE_ARGS NAMES Catch2)
  FetchContent_MakeAvailable(Catch2)
  include(Catch)
endif()

find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
  message(STATUS "Python3 found: ${Python3_EXECUTABLE}")
  # Probe for third-party packages used by cincuenta test scripts.
  # Add one block per package; gate the relevant tests on the result.
  foreach(_pkg IN ITEMS numpy)
    execute_process(
      COMMAND ${Python3_EXECUTABLE} -c "import ${_pkg}"
      RESULT_VARIABLE _rc
      OUTPUT_QUIET ERROR_QUIET)
    if(_rc EQUAL 0)
      set(Python3_${_pkg}_FOUND TRUE)
      message(STATUS "  Python3 package ${_pkg}: found")
    else()
      set(Python3_${_pkg}_FOUND FALSE)
      message(STATUS "  Python3 package ${_pkg}: NOT found — dependent tests skipped")
    endif()
  endforeach()
else()
  message(STATUS "Python3 not found — Python-based tests will be skipped")
endif()

include(GNUInstallDirs)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(APPLE)
  set(CMAKE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
else()
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()

option(DMRG_BUILD_BATCHED "Enable building with Batched support" ON)

include(CMakeDependentOption)
cmake_dependent_option(DMRG_BUILD_BATCHED_MAGMA "Enable building with Magma" OFF "DMRG_BUILD_BATCHED" OFF)

add_subdirectory(util)

add_subdirectory(PsimagLite)
add_subdirectory(dmrg)

add_subdirectory(cincuenta/src)
add_subdirectory(LanczosPlusPlus/src)

option(BUILD_DOC "Build PDF manual (requires pdflatex, bibtex, and perl)" OFF)
if(BUILD_DOC)
  add_subdirectory(doc)
endif()
