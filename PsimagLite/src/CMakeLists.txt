# Add -Wextra in the future
add_library(
  psimaglite
  MersenneTwister.cpp
  Matrix.cpp
  Mpi.cpp
  Concurrency.cpp
  ProgressIndicator.cpp
  MemResolv.cpp
  PsimagLite.cpp
  PsiBase64.cpp
  RedirectOutput.cpp
  SpecialFunctions.cpp
  Io/TypeToH5.cpp
  TridiagonalMatrix.cpp
  PredicateSimple.cpp
  Ainur/AinurSpirit.cpp
  Ainur/AinurConvert.cpp)
target_include_directories(psimaglite PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/Io
                                             ${CMAKE_CURRENT_SOURCE_DIR}/Ainur)

target_compile_features(psimaglite PUBLIC cxx_std_17)
target_link_libraries(psimaglite PUBLIC loki)

add_library(psimaglite::psimaglite ALIAS psimaglite)

# find_package(MPI COMPONENTS CXX)
if(MPI_FOUND)
  target_compile_definitions(psimaglite PUBLIC USE_MPI)
  target_link_libraries(psimaglite PUBLIC MPI::MPI_CXX)

  # Prefer HDF5 parallel when MPI is enabled
  set(HDF5_PREFER_PARALLEL ON)
endif()

find_package(
  HDF5 REQUIRED
  COMPONENTS CXX HL
  MODULE)
if(MPI_FOUND AND NOT HDF5_IS_PARALLEL)
  message(WARNING "HDF5 that was found \(${HDF5_DIR}\) is not parallel")
endif()
target_link_libraries(psimaglite PUBLIC HDF5::HDF5)

find_package(Boost)
if(Boost_FOUND)
  target_compile_definitions(psimaglite PUBLIC USE_BOOST)
  target_link_libraries(psimaglite PUBLIC Boost::boost)
endif()

find_package(Threads)
if(Threads_FOUND)
  target_compile_definitions(psimaglite PUBLIC USE_PTHREADS)
  target_link_libraries(psimaglite PUBLIC Threads::Threads)
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(psimaglite PUBLIC BLAS::BLAS)
target_link_libraries(psimaglite PUBLIC LAPACK::LAPACK)

find_package(GSL REQUIRED)
target_compile_definitions(psimaglite PUBLIC USE_GSL)
target_link_libraries(psimaglite PUBLIC GSL::gsl)

if(BUILD_SHARED_LIBS)
  install(TARGETS psimaglite LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

if(BUILD_TESTING)
  add_subdirectory(tests)

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Ainur/function.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
  add_executable(Ainur_Test1 Ainur/test.cpp)
  target_link_libraries(Ainur_Test1 psimaglite::psimaglite)
  add_test(NAME Ainur_Test1 COMMAND Ainur_Test1 ${CMAKE_CURRENT_SOURCE_DIR}/Ainur/prelude.ain
                                    ${CMAKE_CURRENT_SOURCE_DIR}/Ainur/input0.ain)

  add_executable(Ainur_Test2 Ainur/testReal.cpp)
  target_link_libraries(Ainur_Test2 psimaglite::psimaglite)
  add_test(NAME Ainur_Test2 COMMAND Ainur_Test2 ${CMAKE_CURRENT_SOURCE_DIR}/Ainur/input1.ain)
endif()
