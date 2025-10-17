
add_library(
  psimaglite
  STATIC
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
  psimaglite PUBLIC ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/Io
                    ${CMAKE_SOURCE_DIR}/Ainur)
target_link_libraries(psimaglite PRIVATE HDF5::HDF5)
