set(SRC_FILES
    apply_Htarget_vbatch.cpp
    cal_kron_flops.cpp
    dmrg_lapack.cpp
    dmrg_malloc.cpp
    dmrg_vbatch.cpp
    gen_patches_comb.cpp
    setup_matrix.cpp
    setup_vbatch.cpp
    unsetup_vbatch.cpp
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
endif()

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
