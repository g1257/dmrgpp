# Add -Wextra in the future
add_library(
  psimaglite
  PsimagLite/MersenneTwister.cpp
  PsimagLite/Matrix.cpp
  PsimagLite/Mpi.cpp
  PsimagLite/Concurrency.cpp
  PsimagLite/ProgressIndicator.cpp
  PsimagLite/MemResolv.cpp
  PsimagLite/PsimagLite.cpp
  PsimagLite/PsiBase64.cpp
  PsimagLite/SpecialFunctions.cpp
  PsimagLite/Io/TypeToH5.cpp
  PsimagLite/TridiagonalMatrix.cpp
  PsimagLite/PredicateSimple.cpp
  PsimagLite/Ainur/AinurSpirit.cpp
  PsimagLite/Ainur/AinurConvert.cpp)

# Resolve Kokkos and KokkosKernels as a coherent pair. An installed
# KokkosKernels package must be considered before a Kokkos target exists so its
# package configuration can select the Kokkos installation it was built with.
include(FetchContent)

if(TARGET Kokkos::kokkoskernels)
  message(
    FATAL_ERROR
      "A pre-existing Kokkos::kokkoskernels target cannot be verified as coherent with Kokkos::kokkos; let DMRG++ resolve the pair"
  )
endif()

set(_dmrgpp_probe_kokkoskernels TRUE)
if(FETCHCONTENT_TRY_FIND_PACKAGE_MODE STREQUAL "NEVER"
   OR FETCHCONTENT_SOURCE_DIR_KOKKOS
   OR FETCHCONTENT_SOURCE_DIR_KOKKOSKERNELS
   OR TARGET Kokkos::kokkos)
  set(_dmrgpp_probe_kokkoskernels FALSE)
endif()

if(_dmrgpp_probe_kokkoskernels)
  list(APPEND CMAKE_MESSAGE_INDENT "[KokkosKernels] ")
  find_package(KokkosKernels 5.1 QUIET CONFIG NAMES KokkosKernels)
  list(POP_BACK CMAKE_MESSAGE_INDENT)
endif()

if(KokkosKernels_FOUND)
  if(NOT TARGET Kokkos::kokkos OR NOT TARGET Kokkos::kokkoskernels)
    message(FATAL_ERROR "The KokkosKernels package did not provide both Kokkos::kokkos and Kokkos::kokkoskernels")
  endif()
  set(_dmrgpp_kokkos_provenance "package via KokkosKernels")
  set(_dmrgpp_kokkoskernels_provenance "package")

  # Remove values forced for a source build by an earlier configure. Newer
  # installed packages export these as normal variables that can be checked.
  unset(KokkosKernels_ENABLE_TPL_BLAS CACHE)
  unset(KokkosKernels_ENABLE_TPL_LAPACK CACHE)
  if(KokkosKernels_VERSION VERSION_GREATER_EQUAL 5.2.99)
    if(NOT KokkosKernels_ENABLE_TPL_BLAS)
      message(FATAL_ERROR "KokkosKernels must be configured with BLAS support!")
    endif()
    if(NOT KokkosKernels_ENABLE_TPL_LAPACK)
      message(FATAL_ERROR "KokkosKernels must be configured with LAPACK support!")
    endif()
  endif()
else()
  # No installed pair was selected. Resolve Kokkos on its own, then always
  # build KokkosKernels against that selected target. Package discovery for
  # KokkosKernels is disabled so an installed library cannot be substituted.
  if(TARGET Kokkos::kokkos)
    set(_dmrgpp_kokkos_provenance "provided")
  else()
    list(APPEND CMAKE_MESSAGE_INDENT "[Kokkos] ")
    FetchContent_Declare(
      Kokkos
      GIT_REPOSITORY https://github.com/kokkos/kokkos
      GIT_TAG 5.2.0
      SYSTEM FIND_PACKAGE_ARGS 5.1 NAMES Kokkos)
    FetchContent_MakeAvailable(Kokkos)
    FetchContent_GetProperties(Kokkos SOURCE_DIR _dmrgpp_kokkos_source_dir)
    list(POP_BACK CMAKE_MESSAGE_INDENT)

    if(NOT TARGET Kokkos::kokkos)
      message(FATAL_ERROR "Kokkos did not provide the required Kokkos::kokkos target")
    endif()
    if(_dmrgpp_kokkos_source_dir)
      set(_dmrgpp_kokkos_provenance "source")
    else()
      set(_dmrgpp_kokkos_provenance "package or provider")
    endif()
  endif()

  # Since DMRG++ requires BLAS and LAPACK, require the source build to use both
  # as TPLs.
  set(KokkosKernels_ENABLE_TPL_BLAS
      ON
      CACHE BOOL "Whether to enable the BLAS TPL" FORCE)
  set(KokkosKernels_ENABLE_TPL_LAPACK
      ON
      CACHE BOOL "Whether to enable the LAPACK TPL" FORCE)
  list(APPEND CMAKE_MESSAGE_INDENT "[KokkosKernels] ")
  block()
  set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE NEVER)
  FetchContent_Declare(
    KokkosKernels
    GIT_REPOSITORY https://github.com/kokkos/kokkos-kernels
    GIT_TAG 5.2.0
    SYSTEM)
  endblock()
  FetchContent_MakeAvailable(KokkosKernels)
  FetchContent_GetProperties(KokkosKernels SOURCE_DIR _dmrgpp_kokkoskernels_source_dir)
  list(POP_BACK CMAKE_MESSAGE_INDENT)

  if(NOT TARGET Kokkos::kokkoskernels)
    message(FATAL_ERROR "KokkosKernels did not provide the required Kokkos::kokkoskernels target")
  endif()
  if(_dmrgpp_kokkoskernels_source_dir)
    set(_dmrgpp_kokkoskernels_provenance "source")
  else()
    set(_dmrgpp_kokkoskernels_provenance "provider")
  endif()
endif()

set(DMRGPP_KOKKOS_PROVENANCE
    "${_dmrgpp_kokkos_provenance}"
    CACHE INTERNAL "How Kokkos was resolved" FORCE)
set(DMRGPP_KOKKOSKERNELS_PROVENANCE
    "${_dmrgpp_kokkoskernels_provenance}"
    CACHE INTERNAL "How KokkosKernels was resolved" FORCE)
message(
  STATUS "Kokkos dependency pair: Kokkos=${DMRGPP_KOKKOS_PROVENANCE}, KokkosKernels=${DMRGPP_KOKKOSKERNELS_PROVENANCE}")

target_link_libraries(psimaglite PUBLIC Kokkos::kokkoskernels Kokkos::kokkos)
target_sources(psimaglite PRIVATE PsimagLite/KokkosGemm.cpp)

target_include_directories(psimaglite PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/Io
                                             ${CMAKE_CURRENT_SOURCE_DIR}/Ainur)

target_compile_features(psimaglite PUBLIC cxx_std_17)
target_link_libraries(psimaglite PUBLIC loki)

add_library(psimaglite::psimaglite ALIAS psimaglite)

if(DMRGPP_USE_MPI)
  target_compile_definitions(psimaglite PUBLIC DMRGPP_USE_MPI)
  target_link_libraries(psimaglite PUBLIC MPI::MPI_CXX)
endif()

# Prefer HDF5 parallel when MPI is enabled.
set(HDF5_PREFER_PARALLEL ${DMRGPP_USE_MPI})
find_package(
  HDF5 REQUIRED
  COMPONENTS CXX HL
  MODULE)
if(DMRGPP_USE_MPI AND NOT HDF5_IS_PARALLEL)
  message(WARNING "HDF5 that was found \(${HDF5_DIR}\) is not parallel")
endif()
target_link_libraries(psimaglite PUBLIC HDF5::HDF5)

find_package(Boost REQUIRED)
target_link_libraries(psimaglite PUBLIC Boost::boost)

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
target_link_libraries(psimaglite PUBLIC GSL::gsl)

if(BUILD_SHARED_LIBS)
  install(TARGETS psimaglite LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

if(BUILD_TESTING)
  add_subdirectory(tests)

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/PsimagLite/Ainur/function.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
  add_executable(Ainur_Test1 PsimagLite/Ainur/test.cpp)
  target_link_libraries(Ainur_Test1 psimaglite::psimaglite)
  add_test(NAME Ainur_Test1 COMMAND Ainur_Test1 ${CMAKE_CURRENT_SOURCE_DIR}/PsimagLite/Ainur/prelude.ain
                                    ${CMAKE_CURRENT_SOURCE_DIR}/PsimagLite/Ainur/input0.ain)

  add_executable(Ainur_Test2 PsimagLite/Ainur/testReal.cpp)
  target_link_libraries(Ainur_Test2 psimaglite::psimaglite)
  add_test(NAME Ainur_Test2 COMMAND Ainur_Test2 ${CMAKE_CURRENT_SOURCE_DIR}/PsimagLite/Ainur/input1.ain)
endif()
