add_executable(test_KadanoffBaym test_KadanoffBaym.cpp)
target_link_libraries(test_KadanoffBaym PRIVATE Catch2::Catch2WithMain psimaglite)
target_include_directories(test_KadanoffBaym PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src)
catch_discover_tests(test_KadanoffBaym)

add_executable(test_spectralWeight test_spectralWeight.cpp)
target_link_libraries(test_spectralWeight PRIVATE Catch2::Catch2WithMain psimaglite dmrgpp_utils)
target_include_directories(test_spectralWeight PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src
                                                       ${CMAKE_SOURCE_DIR}/dmrg/Engine)
catch_discover_tests(test_spectralWeight)

add_executable(test_ImpuritySolverNeqExactDiag test_ImpuritySolverNeqExactDiag.cpp)
target_link_libraries(test_ImpuritySolverNeqExactDiag PRIVATE Catch2::Catch2WithMain dmrgpp_utils lanczos_globals
                                                              dmrg_runner kronutil)
target_include_directories(
  test_ImpuritySolverNeqExactDiag PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src ${CMAKE_SOURCE_DIR}/dmrg/Engine
                                          ${CMAKE_SOURCE_DIR}/LanczosPlusPlus/src/Engine ${CMAKE_SOURCE_DIR})
catch_discover_tests(test_ImpuritySolverNeqExactDiag)
