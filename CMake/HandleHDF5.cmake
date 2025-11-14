# Copyright (c) 2025 , UT-Battelle, LLC All rights reserved

# UT Battelle Open Source Software License 11242008 See LICENSE.txt

# Parts of this file are adapted from QMCPACK
# http://https://github.com/QMCPACK/qmcpack

# Licensed under the University of Illinois/NCSA Open Source License.

if(MPI_FOUND)
  option(HDF5_PREFER_PARALLEL "Request parallel/serial HDF5 library" ON)
else(MPI_FOUND)
  option(HDF5_PREFER_PARALLEL "Request parallel/serial HDF5 library" OFF)
  if(HDF5_PREFER_PARALLEL)
    message(FATAL_ERROR "Parallel HDF5 library cannot be selected with QMCPACK non-MPI build. "
                        "Please set HDF5_PREFER_PARALLEL=0.")
  endif(HDF5_PREFER_PARALLEL)
endif(MPI_FOUND)

if(PSIMAGLITE_BUILD_STATIC)
  message(STATUS "Linking static HDF5 library")
  set(HDF5_USE_STATIC_LIBRARIES on)
else()
  message(STATUS "Linking dynamic HDF5 library")
  set(HDF5_USE_STATIC_LIBRARIES off)
endif()

find_package(
  HDF5
  COMPONENTS C CXX HL
  MODULE) # Note: minimum version check is done
# below to bypass find_package
# and HDF5 version compatibility subtleties

if(HDF5_FOUND)
  if(HDF5_VERSION)
    if(HDF5_VERSION VERSION_LESS 1.10.0)
      message(FATAL_ERROR "PsimagLite requires HDF5 version >= 1.10.0")
    endif()
  endif(HDF5_VERSION)
  if(HDF5_IS_PARALLEL)
    if(MPI_FOUND)
      message(STATUS "Parallel HDF5 library found")
      option(ENABLE_PHDF5 "Enable code paths using parallel HDF5" ON)
    else(MPI_FOUND)
      message(FATAL_ERROR "Parallel HDF5 library found but cannot be used with PsimagLite non-MPI build. "
                          "Please provide a serial HDF5 library or switch to building PsimagLite with MPI.")
    endif(MPI_FOUND)
  else(HDF5_IS_PARALLEL)
    message(STATUS "Serial HDF5 library found")
    option(ENABLE_PHDF5 "Enable code paths using parallel HDF5" OFF)
    if(ENABLE_PHDF5)
      if(MPI_FOUND)
        message(FATAL_ERROR "Parallel HDF5 code paths requested but serial HDF5 library found! "
                            "Please either provide parallel HDF5 library or set ENABLE_PHDF5=0.")
      else(MPI_FOUND)
        message(FATAL_ERROR "Parallel HDF5 code paths cannot be enabled on non-MPI builds! Please set ENABLE_PHDF5=0.")
      endif(MPI_FOUND)
    endif(ENABLE_PHDF5)
  endif(HDF5_IS_PARALLEL)

  if(ENABLE_PHDF5)
    message(STATUS "Using HDF5 parallel collective I/O code paths")
  else(ENABLE_PHDF5)
    message(STATUS "Using HDF5 non-scalable serial I/O code paths")
  endif(ENABLE_PHDF5)

  if(MPI_FOUND AND NOT ENABLE_PHDF5)
    message(
      WARNING
        "MPI builds may have performance loss by not using parallel HDF5! (Safe to ignore for workstation builds).")
  endif()

  if(CMAKE_BUILD_TYPE AND HDF5_LIBRARIES_DEBUG)
    if(CMAKE_BUILD_TYPE MATCHES DEBUG)
      set(HDF5_LIBRARIES ${HDF5_LIBRARIES_DEBUG})
    else()
      set(HDF5_LIBRARIES ${HDF5_LIBRARIES_RELEASE})
    endif()
  endif()

else(HDF5_FOUND)
  message(FATAL_ERROR "HDF5 not found. Set HDF5_ROOT")
endif(HDF5_FOUND)
