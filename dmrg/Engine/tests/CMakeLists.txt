add_executable(test_LastKrylovSlots test_LastKrylovSlots.cpp)
target_link_libraries(test_LastKrylovSlots PRIVATE Catch2::Catch2WithMain psimaglite)
target_include_directories(test_LastKrylovSlots PRIVATE ${CMAKE_SOURCE_DIR}/dmrg/Engine)
catch_discover_tests(test_LastKrylovSlots)

# BatchedGemm test
if(DMRG_BUILD_BATCHED)
  add_executable(test_BatchedGemm test_BatchedGemm.cpp)
  target_link_libraries(test_BatchedGemm PRIVATE Catch2::Catch2WithMain dmrgpp_utils psimaglite gpuplugin)
  target_include_directories(test_BatchedGemm PRIVATE ${CMAKE_SOURCE_DIR}/dmrg/Engine
                                                      ${CMAKE_SOURCE_DIR}/dmrg/Engine/MatrixVectorKron)
  catch_discover_tests(test_BatchedGemm)
endif()
