add_executable(test_LastKrylovSlots test_LastKrylovSlots.cpp)
target_link_libraries(test_LastKrylovSlots PRIVATE Catch2::Catch2WithMain psimaglite)
target_include_directories(test_LastKrylovSlots PRIVATE ${CMAKE_SOURCE_DIR}/dmrg/Engine)
catch_discover_tests(test_LastKrylovSlots)

add_executable(testVectorWithOffsets testVectorWithOffsets.cpp)
target_link_libraries(testVectorWithOffsets PRIVATE Catch2::Catch2WithMain dmrgpp_utils psimaglite)
target_include_directories(testVectorWithOffsets PRIVATE ${CMAKE_SOURCE_DIR}/dmrg ${CMAKE_SOURCE_DIR}/dmrg/Engine)
catch_discover_tests(testVectorWithOffsets)

# BatchedGemm test
add_executable(test_BatchedGemm test_BatchedGemm.cpp)
target_link_libraries(test_BatchedGemm PRIVATE Catch2::Catch2WithMain dmrgpp_utils psimaglite)
if(PLUGIN_SC)
  target_link_libraries(test_BatchedGemm PRIVATE gpuplugin)
endif()
target_include_directories(test_BatchedGemm PRIVATE ${CMAKE_SOURCE_DIR}/dmrg/Engine
                                                    ${CMAKE_SOURCE_DIR}/dmrg/Engine/MatrixVectorKron)
catch_discover_tests(test_BatchedGemm)

add_executable(test_LdaggerL test_LdaggerL.cpp)
target_link_libraries(test_LdaggerL PRIVATE Catch2::Catch2WithMain)
target_compile_definitions(
  test_LdaggerL PRIVATE DMRG_EXECUTABLE_PATH="$<TARGET_FILE:dmrg>"
                        DMRG_LIOUVILLIAN_INPUT_PATH="${CMAKE_SOURCE_DIR}/TestSuite/inputs/input8902.ain")
add_dependencies(test_LdaggerL dmrg)
catch_discover_tests(test_LdaggerL)
