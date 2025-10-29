cmake_minimum_required(VERSION 3.21)

find_package(BLAS)
find_package(LAPACK)

# Not sure where to put the build types
if (NOT  ${CMAKE_BUILD_TYPE}  STREQUAL "")
  string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)
  if (build_type STREQUAL debug)
    add_compile_options(-g3 -Werror)
  else()
    add_compile_options(-O3 -DNDEBUG)
  endif()
endif()

# testLanczos executable
add_executable(testLanczos
               testLanczos.cpp)
target_include_directories(testLanczos
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/../src
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/..)

target_link_libraries(
  testLanczos
  psimaglite
  BLAS::BLAS
  LAPACK::LAPACK
)
