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
  # 1. Load the PkgConfig module
  find_package(PkgConfig REQUIRED)

  set(ENV{PKG_CONFIG_PATH} "/home/svh/Documents/CONNEQT/magma-install/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")

  # 2. Check for the module and create an imported target
  # Example: searching for 'glib-2.0'
  pkg_check_modules(MAGMA REQUIRED IMPORTED_TARGET magma)
endif()

#find_package(OpenMP REQUIRED COMPONENTS C)

# 3. Link against the created target
# The target is automatically named PkgConfig::<PREFIX>
#target_link_libraries(my_app PRIVATE PkgConfig::GLIB2)

#set(magma_DIR "/home/svh/Documents/CONNEQT/magma-2.10.0/build")
#find_package(magma REQUIRED)

add_library(gpuplugin ${SRC_FILES})
target_include_directories(gpuplugin PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(gpuplugin PUBLIC psimaglite::psimaglite PkgConfig::MAGMA)
target_compile_definitions(gpuplugin PUBLIC FpType=double USE_MAGMA)
if(DMRG_BUILD_BATCHED_MAGMA)
  target_compile_definitions(gpuplugin PUBLIC USE_MAGMA)
endif()

add_executable(test_vbatch test_vbatch.cpp)
target_link_libraries(test_vbatch gpuplugin)
add_test(NAME vbatch COMMAND test_vbatch)

add_executable(test_gen_patches_comb test_gen_patches_comb.cpp)
target_link_libraries(test_gen_patches_comb gpuplugin)
add_test(NAME gen_patches_comb COMMAND test_gen_patches_comb)
