set(SRC_FILES
    apply_Htarget_vbatch.cpp
    cal_kron_flops.cpp
    dmrg_malloc.cpp
    dmrg_vbatch.cpp
    gen_patches_comb.cpp
    setup_matrix.cpp
    setup_vbatch.cpp
    unsetup_vbatch.cpp
    setup_nC.cpp
    setup_sparse_batch.cpp
    unsetup_sparse_batch.cpp
    apply_Htarget_sparse.cpp
    estimate_work.cpp
    BatchedGemm.cpp)

add_library(gpuplugin ${SRC_FILES})
target_include_directories(gpuplugin PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}) # ${PROJECT_SOURCE_DIR}/PsimagLite/src)
target_link_libraries(gpuplugin PUBLIC psimaglite::psimaglite) # HDF5::HDF5)
target_compile_definitions(gpuplugin PUBLIC FpType=double)
