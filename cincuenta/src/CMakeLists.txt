# This is repeated here but should possibly just use the custom
# command and info from dmrgpp/src

add_executable(cincuenta cincuenta.cpp)
target_link_libraries(cincuenta PUBLIC dmrgpp_utils lanczos_globals dmrg_runner kronutil)

add_executable(freeGimp freeGimp.cpp)
target_link_libraries(freeGimp PUBLIC dmrgpp_utils kronutil)

install(TARGETS cincuenta freeGimp RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Build the unit tests
# Missing for old work still minimal for new dev
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()

# Executable Level Tests
set(PATH_TO_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite/inputs)

# exact impurity solver
add_test(NAME cincuenta0 COMMAND ./cincuenta -I ${PATH_TO_INPUTS} -f inputCincuenta0.ain) # energy 0.456977
# dmrg impurity solver, 3-bath fast CI check — energy -4.86046
add_test(NAME cincuentaDmrg3Bath COMMAND ./cincuenta -I ${PATH_TO_INPUTS} -f inputCincuentaDmrg3Bath.ain)
# dmrg impurity solver — slow; excluded from fast CI, run with: ctest -L Nightly
add_test(NAME cincuenta1 COMMAND ./cincuenta -I ${PATH_TO_INPUTS} -f inputCincuenta1.ain) # energy 0.456977
set_tests_properties(cincuenta1 PROPERTIES LABELS "Nightly")

# Production GBEK must not silently fall back to the obsolete fixed first-bath
# approximation when NeqBathRank is omitted or explicitly zero.  The wrapper
# treats the program's intentional uncaught error as the expected result.
if(Python3_FOUND)
  add_test(
    NAME neqGBEKRejectRankZero
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite/check_neq_gbek_rank_zero_rejected.py
            $<TARGET_FILE:cincuenta> ${PATH_TO_INPUTS})
  add_test(
    NAME neqGBEKRejectIterZero
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite/check_neq_gbek_iter_zero_rejected.py
            $<TARGET_FILE:cincuenta> ${PATH_TO_INPUTS})
  add_test(
    NAME neqGBEKRejectRemovedSolver
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite/check_neq_gbek_removed_solver_rejected.py
            $<TARGET_FILE:cincuenta> ${PATH_TO_INPUTS})
  add_test(
    NAME neqGBEKRejectRemovedTolerance
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite/check_neq_gbek_removed_tolerance_rejected.py
            $<TARGET_FILE:cincuenta> ${PATH_TO_INPUTS})
endif()

# ---- U=0 self-consistency tests: G_imp must equal G_cluster_0 ----------------
# Physical basis: at U=0 the self-energy Sigma=0 exactly, so the DMFT condition
# G_imp = local lattice GF must hold for any correct impurity solver.
# cincuenta writes gimp_<solver>.txt and latticeG_<solver>.txt after the Matsubara
# loop; compare_gimp_latticeG.py checks they agree to within bath-fitting error.
if(Python3_FOUND)
  if(NOT PATH_TO_TESTSUITE)
    set(PATH_TO_TESTSUITE ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite)
  endif()

  # -- ExactDiag ---------------------------------------------------------------
  set(U0_ED_DIR ${CMAKE_CURRENT_BINARY_DIR}/u0ExactDiag)
  file(MAKE_DIRECTORY ${U0_ED_DIR})

  add_test(
    NAME u0_gImpEqGcluster_exactdiag
    COMMAND $<TARGET_FILE:cincuenta> -I ${PATH_TO_INPUTS} -f inputU0ExactDiag.ain
    WORKING_DIRECTORY ${U0_ED_DIR})

  add_test(NAME u0_gImpEqGcluster_exactdiag_check
           COMMAND ${Python3_EXECUTABLE} ${PATH_TO_TESTSUITE}/compare_gimp_latticeG.py ${U0_ED_DIR}/gimp_exactdiag.txt
                   ${U0_ED_DIR}/latticeG_exactdiag.txt --tolerance 0.015)
  set_tests_properties(u0_gImpEqGcluster_exactdiag_check PROPERTIES DEPENDS u0_gImpEqGcluster_exactdiag)

  # -- DMRG --------------------------------------------------------------------
  set(U0_DMRG_DIR ${CMAKE_CURRENT_BINARY_DIR}/u0Dmrg)
  file(MAKE_DIRECTORY ${U0_DMRG_DIR})

  add_test(
    NAME u0_gImpEqGcluster_dmrg
    COMMAND $<TARGET_FILE:cincuenta> -I ${PATH_TO_INPUTS} -f inputU0Dmrg.ain
    WORKING_DIRECTORY ${U0_DMRG_DIR})

  add_test(NAME u0_gImpEqGcluster_dmrg_check
           COMMAND ${Python3_EXECUTABLE} ${PATH_TO_TESTSUITE}/compare_gimp_latticeG.py ${U0_DMRG_DIR}/gimp_dmrg.txt
                   ${U0_DMRG_DIR}/latticeG_dmrg.txt --tolerance 0.05)
  set_tests_properties(u0_gImpEqGcluster_dmrg_check PROPERTIES DEPENDS u0_gImpEqGcluster_dmrg LABELS "Nightly")
  set_tests_properties(u0_gImpEqGcluster_dmrg PROPERTIES LABELS "Nightly")

endif()

# Fast end-to-end GBEK regression: a fitted non-atomic first bath plus a
# positive-rank (L=2) Cholesky second bath.  Keep the output isolated so this
# test cannot overwrite other non-equilibrium runs.
if(Python3_FOUND)
  set(NEQ_GBEK_FAST_DIR ${CMAKE_CURRENT_BINARY_DIR}/neqGbekFast)
  file(MAKE_DIRECTORY ${NEQ_GBEK_FAST_DIR})

  add_test(
    NAME neqGBEKFast
    COMMAND $<TARGET_FILE:cincuenta> -I ${PATH_TO_INPUTS} -f inputNeqGBEKFast.ain
    WORKING_DIRECTORY ${NEQ_GBEK_FAST_DIR})
  set_tests_properties(neqGBEKFast PROPERTIES FIXTURES_SETUP neq_gbek_fast)

  add_test(
    NAME neqGBEKFast_check
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite/check_neq_gbek_fast.py
            ${NEQ_GBEK_FAST_DIR} --matsubaras 16)
  set_tests_properties(neqGBEKFast_check PROPERTIES FIXTURES_REQUIRED neq_gbek_fast)
endif()
