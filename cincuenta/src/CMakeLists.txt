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

# General helper: check one captured regex group against an expected value.
# check_value uses the LAST match by default; pass --first as part of pattern
# arg (prepend to pattern string) to stop at the first match instead.
function(
  add_check
  name
  pattern
  expected
  tolerance
  file
  parent_test)
  add_test(NAME ${name} COMMAND $<TARGET_FILE:check_value> "${pattern}" "${expected}" "${tolerance}" "${file}")
  set_tests_properties(${name} PROPERTIES DEPENDS ${parent_test})
endfunction()

# exact impurity solver
add_test(NAME cincuenta0 COMMAND ./cincuenta -I ${PATH_TO_INPUTS} -f inputCincuenta0.ain) # energy 0.456977
# dmrg impurity solver, 3-bath fast CI check — energy -4.86046
add_test(NAME cincuentaDmrg3Bath COMMAND ./cincuenta -I ${PATH_TO_INPUTS} -f inputCincuentaDmrg3Bath.ain)
# dmrg impurity solver — slow; excluded from fast CI, run with: ctest -L Nightly
add_test(NAME cincuenta1 COMMAND ./cincuenta -I ${PATH_TO_INPUTS} -f inputCincuenta1.ain) # energy 0.456977
set_tests_properties(cincuenta1 PROPERTIES LABELS "Nightly")

# Non-equilibrium DMFT: interaction quench U_i=0 -> U_f=2, 5-bath exact diag,
# beta=16, t_max=5 (N_t=200).  Checking correct initial values the
# Bethe lattice; plus a regression anchor for the 5-bath finite-bath discretisation.
add_test(NAME neqTsuji0 COMMAND $<TARGET_FILE:cincuenta> -I ${PATH_TO_INPUTS} -f inputNeqTsuji.ain)
set_tests_properties(neqTsuji0 PROPERTIES LABELS "Nightly")

# G^R(t=0, t'=0): boundary condition G^R(t,t) = -i
add_check(neqTsuji0_GR00_Im "0\\.0000000000 0\\.0000000000 [^ ]+ ([-0-9.eE+]+)" "-1.0" "1e-8"
          ${CMAKE_CURRENT_BINARY_DIR}/green-retarded neqTsuji0)
set_tests_properties(neqTsuji0_GR00_Im PROPERTIES LABELS "Nightly")

# G^M(tau=0) = -n_imp: PH-symmetric bath gives exact half-filling
add_check(neqTsuji0_GM0_Re "^0\\.0000000000 ([-0-9.eE+]+) [-0-9.eE+]+" "-0.5" "1e-8"
          ${CMAKE_CURRENT_BINARY_DIR}/green-matsubara-t neqTsuji0)
set_tests_properties(neqTsuji0_GM0_Re PROPERTIES LABELS "Nightly")

# G^<(t=0, t'=0): initial occupation Im[G^<(0,0)] = n(0) = 0.5 for U_i=0 PH-symmetric
add_check(neqTsuji0_GL00_Im "0\\.0000000000 0\\.0000000000 [^ ]+ ([-0-9.eE+]+)" "0.5" "1e-8"
          ${CMAKE_CURRENT_BINARY_DIR}/green-lesser neqTsuji0)
set_tests_properties(neqTsuji0_GL00_Im PROPERTIES LABELS "Nightly")

# G^R(t=1, t'=0): retarded GF one time unit after quench — regression anchor
add_check(neqTsuji0_GR10_Im "1\\.0000000000 0\\.0000000000 [^ ]+ ([-0-9.eE+]+)" "-0.377208" "1e-4"
          ${CMAKE_CURRENT_BINARY_DIR}/green-retarded neqTsuji0)
set_tests_properties(neqTsuji0_GR10_Im PROPERTIES LABELS "Nightly")

# n(t=1) = Im[G^<(1,1)]: PH-symmetric system locks n=0.5 — regression anchor
add_check(neqTsuji0_n1 "1\\.0000000000 1\\.0000000000 [^ ]+ ([-0-9.eE+]+)" "0.5" "1e-8"
          ${CMAKE_CURRENT_BINARY_DIR}/green-lesser neqTsuji0)
set_tests_properties(neqTsuji0_n1 PROPERTIES LABELS "Nightly")

# n(t=5) = Im[G^<(5,5)]: PH-symmetric system locks n=0.5 — regression anchor
add_check(neqTsuji0_n5 "5\\.0000000000 5\\.0000000000 [^ ]+ ([-0-9.eE+]+)" "0.5" "1e-8"
          ${CMAKE_CURRENT_BINARY_DIR}/green-lesser neqTsuji0)
set_tests_properties(neqTsuji0_n5 PROPERTIES LABELS "Nightly")

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

  # ---- U=0 no-quench neq tests: G_imp must equal Weiss field G_0 ---------------
  # Physical basis: at U_i=U_f=0 (no quench) Sigma_imp=0, so Dyson gives G_imp=G_0.
  # cincuenta writes green-{retarded,lesser} (G_imp) and weiss-green-* (G_0);
  # compare_gimp_weiss.py checks they agree.

  # -- Neq ExactDiag -----------------------------------------------------------
  set(U0_NEQ_ED_DIR ${CMAKE_CURRENT_BINARY_DIR}/u0NeqExactDiag)
  file(MAKE_DIRECTORY ${U0_NEQ_ED_DIR})

  add_test(
    NAME u0_neq_gImpEqG0_exactdiag
    COMMAND $<TARGET_FILE:cincuenta> -I ${PATH_TO_INPUTS} -f inputU0NeqExactDiag.ain
    WORKING_DIRECTORY ${U0_NEQ_ED_DIR})

  # But it doesn't really!
  # Tolerance 0.30: the finite discrete bath cannot track the growing entanglement
  # of the impurity with the environment during time evolution — the bath state drifts
  # away from the correct hybridization, an effect identified by Gramsch, Balzer,
  # Eckstein, Kollar (Phys. Rev. B 88, 235106, 2013).  G^R error (~0.17) accumulates
  # over TmaxNeq=2; G^< error (~0.29) additionally reflects the particle-only Matsubara
  # GF in buildMatsubara (hole sector missing).  All known code bugs are fixed; residual
  # error is physics — non-SCF finite-bath ED without GEBK correction cannot do better.
  add_test(
    NAME u0_neq_gImpEqG0_exactdiag_check
    COMMAND
      ${Python3_EXECUTABLE} ${PATH_TO_TESTSUITE}/compare_gimp_weiss.py --retarded-gimp ${U0_NEQ_ED_DIR}/green-retarded
      --retarded-weiss ${U0_NEQ_ED_DIR}/weiss-green-retarded --lesser-gimp ${U0_NEQ_ED_DIR}/green-lesser --lesser-weiss
      ${U0_NEQ_ED_DIR}/weiss-green-lesser --tolerance 0.30)
  set_tests_properties(u0_neq_gImpEqG0_exactdiag_check PROPERTIES DEPENDS u0_neq_gImpEqG0_exactdiag LABELS "Nightly")
  set_tests_properties(u0_neq_gImpEqG0_exactdiag PROPERTIES LABELS "Nightly")

endif()
