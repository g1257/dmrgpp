set(SRC_FILES
    apply_Htarget_vbatch.cpp
    cal_kron_flops.cpp
    dmrg_lapack.cpp
    dmrg_vbatch.cpp
    gen_patches_comb.cpp
    setup_matrix.cpp
    setup_vbatch.cpp
    setup_nC.cpp
    setup_sparse_batch.cpp
    unsetup_sparse_batch.cpp
    apply_Htarget_sparse_not_ready.cpp
    estimate_work.cpp
    BatchedGemm.cpp)

if(DMRG_BUILD_BATCHED_MAGMA)
  #export PKG_CONFIG_PATH=<MAGMA_INSTALL_PREFIX>/lib/pkgconfig:$PKG_CONFIG_PATH

  find_package(PkgConfig REQUIRED)

  pkg_check_modules(MAGMA REQUIRED IMPORTED_TARGET magma)

  if(Kokkos_ENABLE_HIP)
    list(APPEND SRC_FILES "dmrg_malloc_hip.cpp")
  elseif(Kokkos_ENABLE_CUDA)
    list(APPEND SRC_FILES "dmrg_malloc_cuda.cpp")
  else()
    message(FATAL_ERROR "MAGMA requires either the HIP or CUDA backend!")
  endif()
else()
  list(APPEND SRC_FILES "dmrg_malloc.cpp")
endif()

set(FpType "double")
if(DMRG_BUILD_BATCHED_MAGMA)
  set(USE_MAGMA ON)
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/GPUPluginConfig.h.in ${CMAKE_CURRENT_BINARY_DIR}/GPUPluginConfig.h @ONLY)

add_library(gpuplugin ${SRC_FILES})
target_include_directories(gpuplugin PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(gpuplugin PUBLIC psimaglite::psimaglite)
if(DMRG_BUILD_BATCHED_MAGMA)
  target_link_libraries(gpuplugin PUBLIC PkgConfig::MAGMA)
endif()

add_executable(test_vbatch test_vbatch.cpp)
target_link_libraries(test_vbatch gpuplugin)
add_test(NAME vbatch COMMAND test_vbatch)

add_executable(test_gen_patches_comb test_gen_patches_comb.cpp)
target_link_libraries(test_gen_patches_comb gpuplugin)
add_test(NAME gen_patches_comb COMMAND test_gen_patches_comb)
