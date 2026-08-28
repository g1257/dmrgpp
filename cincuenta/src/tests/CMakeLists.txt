add_executable(test_NeqBathDecomposition test_NeqBathDecomposition.cpp)
target_link_libraries(test_NeqBathDecomposition PRIVATE Catch2::Catch2WithMain
                                                        psimaglite)
target_include_directories(test_NeqBathDecomposition
                           PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src)
catch_discover_tests(test_NeqBathDecomposition)

add_executable(test_NeqRealTimeGf test_NeqRealTimeGf.cpp)
target_link_libraries(test_NeqRealTimeGf PRIVATE Catch2::Catch2WithMain
                                                 psimaglite)
target_include_directories(test_NeqRealTimeGf
                           PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src)
catch_discover_tests(test_NeqRealTimeGf)

add_executable(test_spectralWeight test_spectralWeight.cpp)
target_link_libraries(test_spectralWeight PRIVATE Catch2::Catch2WithMain
                                                  psimaglite dmrgpp_utils)
target_include_directories(
  test_spectralWeight PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src
                              ${CMAKE_SOURCE_DIR}/dmrg/Engine)
catch_discover_tests(test_spectralWeight)

add_executable(test_ImpuritySolverNeqGBEK test_ImpuritySolverNeqGBEK.cpp)
target_link_libraries(
  test_ImpuritySolverNeqGBEK PRIVATE Catch2::Catch2WithMain dmrgpp_utils
                                     lanczos_globals dmrg_runner kronutil)
target_include_directories(
  test_ImpuritySolverNeqGBEK
  PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src ${CMAKE_SOURCE_DIR}/dmrg/Engine
          ${CMAKE_SOURCE_DIR}/LanczosPlusPlus/src/Engine ${CMAKE_SOURCE_DIR})
catch_discover_tests(test_ImpuritySolverNeqGBEK)

add_executable(test_NeqGbekInput test_NeqGbekInput.cpp)
target_link_libraries(
  test_NeqGbekInput PRIVATE Catch2::Catch2WithMain dmrgpp_utils lanczos_globals
                            dmrg_runner kronutil)
target_include_directories(
  test_NeqGbekInput
  PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src ${CMAKE_SOURCE_DIR}/dmrg/Engine
          ${CMAKE_SOURCE_DIR}/LanczosPlusPlus/src/Engine ${CMAKE_SOURCE_DIR})
catch_discover_tests(test_NeqGbekInput)

add_executable(test_DmftSolverEquilibriumInitialData
               test_DmftSolverEquilibriumInitialData.cpp)
target_link_libraries(
  test_DmftSolverEquilibriumInitialData
  PRIVATE Catch2::Catch2WithMain dmrgpp_utils lanczos_globals dmrg_runner
          kronutil)
target_include_directories(
  test_DmftSolverEquilibriumInitialData
  PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src ${CMAKE_SOURCE_DIR}/dmrg/Engine
          ${CMAKE_SOURCE_DIR}/LanczosPlusPlus/src/Engine ${CMAKE_SOURCE_DIR})
set(DMFT_EQUILIBRIUM_INITIAL_DATA_DIR
    ${CMAKE_CURRENT_BINARY_DIR}/equilibriumInitialData)
file(MAKE_DIRECTORY ${DMFT_EQUILIBRIUM_INITIAL_DATA_DIR})
add_test(
  NAME test_DmftSolverEquilibriumInitialData
  COMMAND test_DmftSolverEquilibriumInitialData
  WORKING_DIRECTORY ${DMFT_EQUILIBRIUM_INITIAL_DATA_DIR})
