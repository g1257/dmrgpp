cmake_minimum_required(VERSION 3.22)

project(
  psimaglite
  VERSION 1.1
  LANGUAGES CXX)

set(PROJECT_CMAKE ${psimaglite_SOURCE_DIR}/CMake)
message("Project CMake: ${PROJECT_CMAKE}")
list(APPEND CMAKE_MODULE_PATH ${PROJECT_CMAKE})

set(QMC_CXX_STANDARD
    11
    CACHE STRING "PsiMaglite C++ language standard")
if(NOT QMC_CXX_STANDARD EQUAL 11)
  message(
    WARNING
      "C++11 is the only language standard officially supported by psimaglite. "
      "Using other versions of the C++ standard is unsupported and done entirely at user's own risk."
  )
endif()

# add the automatically determined parts of the RPATH which point to directories
# outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
include(AddRPATH)

include(HandleHDF5)

add_subdirectory(loki)
add_subdirectory(src)
add_subdirectory(examples)
#add_subdirectory(doc)

export(
  EXPORT PsimagLiteTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/PsimagLiteTargets.cmake"
  NAMESPACE psimaglite::)

# Generate config for build tree
file(
  WRITE "${CMAKE_CURRENT_BINARY_DIR}/PsimagLiteConfig.cmake"
  "
include(\"\${CMAKE_CURRENT_LIST_DIR}/PsimagLiteTargets.cmake\")
")

install(
  EXPORT PsimagLiteTargets
  FILE PsimagLiteTargets.cmake
  NAMESPACE psimaglite::
  DESTINATION lib/cmake/psimaglite)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/CMake/PsimagLiteConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/PsimagLiteConfig.cmake"
  INSTALL_DESTINATION lib/cmake/PsimagLite)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/PsimagLiteConfig.cmake"
        DESTINATION lib/cmake/PsimagLite)
