add_executable(test_LastKrylovSlots test_LastKrylovSlots.cpp)
target_link_libraries(test_LastKrylovSlots PRIVATE Catch2::Catch2WithMain psimaglite)
target_include_directories(test_LastKrylovSlots PRIVATE ${CMAKE_SOURCE_DIR}/dmrg/Engine)
catch_discover_tests(test_LastKrylovSlots)

# BatchedGemmKokkos test
add_executable(test_BatchedGemmKokkos test_BatchedGemmKokkos.cpp)
# Link Kokkos and KokkosKernels since the test initializes Kokkos and uses KokkosBatched
target_link_libraries(test_BatchedGemmKokkos PRIVATE Catch2::Catch2WithMain dmrgpp_utils psimaglite Kokkos::kokkos
                                                     KokkosKernels::kokkoskernels gpuplugin)
# Ensure headers from Engine and the MatrixVectorKron subdir are visible
target_include_directories(test_BatchedGemmKokkos PRIVATE ${CMAKE_SOURCE_DIR}/dmrg/Engine
                                                          ${CMAKE_SOURCE_DIR}/dmrg/Engine/MatrixVectorKron)
catch_discover_tests(test_BatchedGemmKokkos)
