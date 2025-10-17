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
target_link_libraries(psimaglite PRIVATE HDF5::HDF5)

export(
  TARGETS psimaglite
  FILE "${CMAKE_CURRENT_BINARY_DIR}/PsimagLite.cmake"
  NAMESPACE PSIMAGLITE::)

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/PsimagLiteConfig.cmake"
     "include(\"\${CMAKE_CURRENT_LIST_DIR}/PsimagLite.cmake\")")

install(
  TARGETS psimaglite
  EXPORT PsimagLite
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)

install(
  EXPORT PsimagLite
  FILE PsimagLite.cmake
  NAMESPACE PSIMAGLITE::
  DESTINATION lib/cmake/psimaglite)

install(DIRECTORY include/ DESTINATION include)
