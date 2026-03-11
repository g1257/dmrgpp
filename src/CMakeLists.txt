add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/GitRevision.h
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/createGitRevision.pl ${CMAKE_CURRENT_BINARY_DIR}/GitRevision.h
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/createGitRevision.pl
  COMMENT "Generating Git Revision File"
  VERBATIM)

add_subdirectory(KronUtil)

# DMRG++ as a library so cincuenta can use it
add_library(dmrgpp ProgramGlobals.cpp Provenance.cpp Utils.cpp Su2Related.cpp Qn.cpp)
target_sources(dmrgpp PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/GitRevision.h)
target_include_directories(dmrgpp PUBLIC ${CMAKE_CURRENT_BINARY_DIR} Engine)
target_link_libraries(dmrgpp PUBLIC psimaglite::psimaglite)

# generate sources file for explicitly instantiating mainLoop4
set(counter 0)
set(instantiationFileList)
foreach(complexOrReal RealType std::complex<RealType>)
  foreach(solverType LanczosSolver ChebyshevSolver)
    foreach(matrixVectorType MatrixVectorOnTheFly MatrixVectorStored MatrixVectorKron)
      math(EXPR counter "${counter} + 1")
      set(instantiationFile ${CMAKE_CURRENT_BINARY_DIR}/DmrgDriver${counter}.cpp)
      configure_file(cmake/DmrgDriver.cpp.in ${instantiationFile})
      list(APPEND instantiationFileList ${instantiationFile})
    endforeach()
  endforeach()
endforeach()
unset(counter)

# Dmrg executable (main executable)
add_executable(dmrg ${instantiationFileList} dmrg.cpp)
target_link_libraries(dmrg dmrgpp kronutil)
target_include_directories(dmrg PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Observe executable
add_executable(observe observe.cpp ObserveDriver0.cpp ObserveDriver1.cpp ObserveDriver2.cpp)
target_link_libraries(observe dmrgpp kronutil)

# Toolbox executable
add_executable(toolboxdmrg toolboxdmrg.cpp Qn.cpp ProgramGlobals.cpp Provenance.cpp)
target_link_libraries(toolboxdmrg dmrgpp)

# Many omegas executable
add_executable(omegas manyOmegas.cpp)
target_link_libraries(omegas dmrgpp kronutil)

if(BUILD_SHARED_LIBS)
  install(TARGETS dmrgpp LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
install(TARGETS dmrg toolboxdmrg omegas RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Tests below for all executables

# For dmrg main executable
set(PATH_TO_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite/inputs)

function(add_lowest_eigenvalue_check_test name expected output_file)
  set(options SHOULD_FAIL)
  set(one_value_args TOLERANCE)
  set(multivalue_args)
  cmake_parse_arguments(PARSE "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

  # name expected output_file)
  if(${AWK} MATCHES ".+-NOTFOUND")
    message(WARNING "No awk found not adding ${name} test")
  else()
    set(tolerance 1e-5)
    if(PARSE_TOLERANCE)
      set(tolerance ${PARSE_TOLERANCE})
    endif()
    add_test(
      NAME ${name}
      COMMAND
        ${AWK} -f ${dmrgpp_SOURCE_DIR}/scripts/check_value.awk -v
        "pattern=Found lowest eigenvalue=[[:space:]]+([.0-9eE+-]+)" -v "expected=${expected}" -v
        "tolerance=${tolerance}" ${CMAKE_CURRENT_BINARY_DIR}/${output_file})
    if(PARSE_SHOULD_FAIL)
      set_tests_properties(${name} PROPERTIES WILL_FAIL TRUE)
    endif()
  endif()

endfunction()

add_test(NAME input2 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input2.ain)
add_lowest_eigenvalue_check_test(output2 -9.51754 runForinput2.cout)
add_test(NAME input3 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input3.ain) # energy -21.5102
add_lowest_eigenvalue_check_test(output3 -21.5102 runForinput3.cout)
add_test(NAME input4 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input4.ain) # energy -15.9938
add_lowest_eigenvalue_check_test(output4 -15.9938 runForinput4.cout)
add_test(NAME input5 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input5.ain) # energy -15.9938
add_lowest_eigenvalue_check_test(output5 -15.9938 runForinput5.cout)
add_test(NAME input7 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input7.ain) # energy -15.9939
add_lowest_eigenvalue_check_test(output7 -15.9939 runForinput7.cout)

# input8 is too slow (disabled)
# add_test(NAME input8 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input8.ain)
add_test(NAME input10 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input10.ain) # energy -19.6747
add_lowest_eigenvalue_check_test(output10 -19.6747 runForinput10.cout)

# input11 and input12 need -DUSE_PTHREADS
# add_test(NAME input11 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input11.ain)
# add_test(NAME input12 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input12.ain)
add_test(NAME input20 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input20.ain) # energy -3.37493
add_lowest_eigenvalue_check_test(output20 -3.37493 runForinput20.cout)
add_test(NAME input21 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input21.ain "<gs|sz|P0>") # none
add_test(NAME input22 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input22.ain "<gs|sz|P0>") # none
add_test(NAME input30 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input30.ain) # energy -4.29307
add_lowest_eigenvalue_check_test(output30 -4.29307 runForinput30.cout)
add_test(NAME input32 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input32.ain) # energy -15.9931
add_lowest_eigenvalue_check_test(output32 -15.9931 runForinput32.cout)
add_test(NAME input60 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input60.ain) # energy -21.7473
add_lowest_eigenvalue_check_test(output60 -21.7473 runForinput60.cout)
add_test(NAME input61 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input61.ain) # energy -5.75007
add_lowest_eigenvalue_check_test(output61 -5.75007 runForinput61.cout)
add_test(NAME input62 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input62.ain) # energy -3.73168
add_lowest_eigenvalue_check_test(output62 -3.73168 runForinput62.cout)
add_test(NAME input63 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input63.ain) # none
add_test(NAME input70 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input70.ain) # energy -9.28874
add_lowest_eigenvalue_check_test(output70 -9.28874 runForinput70.cout)
add_test(NAME input71 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input71.ain) # energy -9.74163
add_lowest_eigenvalue_check_test(output71 -9.74163 runForinput71.cout)
add_test(NAME input80 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input80.ain) # energy -7.37027
add_lowest_eigenvalue_check_test(output80 -7.37027 runForinput80.cout)
add_test(NAME input120 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input120.ain) # energy -9.51754
add_lowest_eigenvalue_check_test(output120 -9.51754 runForinput120.cout)
add_test(NAME input121 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input121.ain -p 12 "<gs|P0>") # none
add_test(NAME input122 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input122.ain -p 12 "<gs|P0>,<gs|P1>,<P0|P1>") # none
add_test(NAME input123 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input123.ain -p 12
                               "<gs|P0>,<P0|P1>,<gs|c'|P2>,<gs|c'|P3>,<P2|P3>") # none
add_test(NAME input126 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input126.ain) # energy -19.6759
add_lowest_eigenvalue_check_test(output126 -19.6759 runForinput126.cout)
add_test(NAME input127 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input127.ain -p 12 "<gs|P0>") # none
add_test(NAME input170 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input170.ain) # energy -4.75877
add_lowest_eigenvalue_check_test(output170 -4.75877 runForinput170.cout)
add_test(NAME input172 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input172.ain) # energy -7.0
add_lowest_eigenvalue_check_test(output172 -7.0 runForinput172.cout)
add_test(NAME input250 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input250.ain) # energy -4.39454
add_lowest_eigenvalue_check_test(output250 -4.39454 runForinput250.cout)
add_test(NAME input320 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input320.ain) # energy -18.8982
add_lowest_eigenvalue_check_test(output320 -18.8982 runForinput320.cout)
add_test(NAME input501 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input501.ain) # energy -6.77334
add_lowest_eigenvalue_check_test(output501 -6.77334 runForinput501.cout)
add_test(NAME input550 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input550.ain) # energy -4.66667
add_lowest_eigenvalue_check_test(output550 -4.66667 runForinput550.cout)
add_test(NAME input1300 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1300.ain) # energy -2.32264
add_lowest_eigenvalue_check_test(output1300 -2.32264 runForinput1300.cout)
add_test(NAME input1302 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1302.ain) # energy -1.78118
add_lowest_eigenvalue_check_test(output1302 -1.78118 runForinput1302.cout)
add_test(NAME input1304 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1304.ain) # energy -1.78118
add_lowest_eigenvalue_check_test(output1304 -1.78118 runForinput1304.cout)
add_test(NAME input1306 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1306.ain) # energy -8.0
add_lowest_eigenvalue_check_test(output1306 -8.0 runForinput1306.cout)
add_test(NAME input1600 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1600.ain) # energy -9.51746

add_lowest_eigenvalue_check_test(output1600 -9.51746 runForinput1600.cout TOLERANCE 1e-4)
add_test(NAME input1814 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1814.ain) # energy -6.0
add_lowest_eigenvalue_check_test(output1814 -6.0 runForinput1814.cout)
add_test(NAME input1815 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1815.ain) # none
add_test(NAME input1817 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1817.ain) # none
# input1818 and 1819 --> too slow (disabled)
# add_test(NAME input1818 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1818.ain)
# add_test(NAME input1819 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1819.ain)
add_test(NAME input2500 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input2500.ain) # energy -4.1632
add_lowest_eigenvalue_check_test(output2500 -4.1632 runForinput2500.cout)
add_test(NAME input2510 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input2510.ain) # none
# needs -DUSE_PTHREADS
# add_test(NAME input2700 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input2700.ain)
add_test(NAME input3020 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input3020.ain) # energy -1.44244
add_lowest_eigenvalue_check_test(output3020 -1.44244 runForinput3020.cout)
add_test(NAME input3021 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input3021.ain "<gs|c'|P2>") # energy -1.44244
add_lowest_eigenvalue_check_test(output3021 -1.44244 runForinput3021.cout)
add_test(NAME input3100 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input3100.ain) # energy -13.9973
add_lowest_eigenvalue_check_test(output3100 -13.9973 runForinput3100.cout)
# SDHS will remain disabled for a while
# add_test(NAME input4650 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input4650.ain)
add_test(NAME input4720 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input4720.ain) # energy 34.0013
add_lowest_eigenvalue_check_test(output4720 34.0013 runForinput4720.cout)
add_test(NAME input6020 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input6020.ain) # energy -2.59474
add_lowest_eigenvalue_check_test(output6020 -2.59474 runForinput6020.cout)
add_test(NAME input6700 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input6700.ain) # energy -4.17889
add_lowest_eigenvalue_check_test(output6700 -4.17889 runForinput6700.cout SHOULD_FAIL)
add_test(NAME input8900 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input8900.ain "<gs|sz_p|gs>") # none
add_test(NAME input8901 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input8901.ain "<gs|sz|gs>") # none
# Tests for the observe executable
add_test(NAME observe0002 COMMAND ./observe -l "+r" -f ${PATH_TO_INPUTS}/input2.ain
                                  "<gs|c';c|gs>,<gs|2.0*sz;2.0*sz|gs>,<gs|n;n|gs>")
set_tests_properties(observe0002 PROPERTIES DEPENDS input2)
add_test(NAME observe0002_1 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/correlations2.pl
                                    ${CMAKE_CURRENT_BINARY_DIR}/runForobserveinput2.cout)
set_tests_properties(observe0002_1 PROPERTIES DEPENDS observe0002)

add_test(NAME observe20 COMMAND ./observe -l "+r" -f ${PATH_TO_INPUTS}/input20.ain "<gs|sz;sz|gs>,<gs|splus;sminus|gs>")
set_tests_properties(observe20 PROPERTIES DEPENDS input20)
add_test(NAME observe20_1 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/correlations2.pl
                                  ${CMAKE_CURRENT_BINARY_DIR}/runForobserveinput20.cout)
set_tests_properties(observe20_1 PROPERTIES DEPENDS observe20)
