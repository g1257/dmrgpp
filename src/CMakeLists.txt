# Add -Wextra in the future
add_library(
  psimaglite STATIC
  MersenneTwister.cpp
  Matrix.cpp
  Mpi.cpp
  Concurrency.cpp
  ProgressIndicator.cpp
  MemResolv.cpp
  PsimagLite.cpp
  PsiBase64.cpp
  SpecialFunctions.cpp
  Io/TypeToH5.cpp
  TridiagonalMatrix.cpp
  PredicateSimple.cpp
  Ainur/AinurSpirit.cpp
  Ainur/AinurConvert.cpp)
target_include_directories(
  psimaglite
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Io>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Ainur>
         $<INSTALL_INTERFACE:include>)

target_link_libraries(psimaglite PUBLIC loki)

set(HDF5_PREFER_PARALLEL ON)
find_package(HDF5 COMPONENTS CXX HL MODULE)
if(HDF5_FOUND AND HDF5_IS_PARALLEL)
  message(STATUS "HDF5 enabled")
  target_link_libraries(psimaglite PUBLIC HDF5::HDF5)
else()
  if(HDF5_FOUND)
    message(WARNING "HDF5 that was found \(${HDF5_DIR}\) is not parallel and cannot be used")
  endif()
  target_compile_definitions(psimaglite PUBLIC USE_IO_SIMPLE)
endif()

find_package(Boost)
if(Boost_FOUND)
  target_compile_definitions(psimaglite PUBLIC USE_BOOST)
  target_link_libraries(psimaglite PUBLIC Boost::boost)
endif()


install(
  TARGETS psimaglite
  EXPORT PsimagLiteTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)

export(
  TARGETS psimaglite
  FILE "${CMAKE_CURRENT_BINARY_DIR}/PsimagLiteTargets.cmake"
  NAMESPACE psimaglite::
  APPEND # Important: append to existing export set
)
