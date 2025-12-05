find_package(Boost REQUIRED CONFIG)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/GitRevision.h
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/createGitRevision.pl ${CMAKE_CURRENT_BINARY_DIR}/GitRevision.h
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/createGitRevision.pl
  COMMENT "Generating Git Revision File"
  VERBATIM)

add_subdirectory(KronUtil)

add_library(Common OBJECT ProgramGlobals.cpp Provenance.cpp Utils.cpp)
target_sources(Common PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/GitRevision.h)
target_link_libraries(Common PRIVATE psimaglite loki HDF5::HDF5)
target_include_directories(Common PRIVATE ${CMAKE_CURRENT_BINARY_DIR} Engine)

add_library(Su2Related OBJECT Su2Related.cpp)
target_link_libraries(Su2Related PRIVATE psimaglite HDF5::HDF5)
target_include_directories(Su2Related PRIVATE ${CMAKE_CURRENT_BINARY_DIR} Engine)

execute_process(
  COMMAND perl -I${CMAKE_CURRENT_SOURCE_DIR} countDrivers.pl
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  ERROR_VARIABLE perl_error
  RESULT_VARIABLE perl_result
  OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REGEX MATCH "([0-9]+)" match_result ${perl_error})
set(NUM_DRIVER_FILES ${CMAKE_MATCH_1})

set(driver_list)
# configure time math string(STRIP "${NUM_DRIVER_FILES}" NUM_DRIVER_FILES_CLEAN)
math(EXPR NUM_DRIVER_FILES_MINUS_ONE "${NUM_DRIVER_FILES} - 1")
foreach(i RANGE 0 ${NUM_DRIVER_FILES_MINUS_ONE})
  list(APPEND driver_list "DmrgDriver${i}.cpp;")
endforeach()

set(driver_templates)
foreach(driver ${driver_list})
  list(APPEND driver_templates "${CMAKE_CURRENT_SOURCE_DIR}/${driver};")
endforeach()
message("driver_list: ${driver_list}")

add_custom_command(
  OUTPUT ${driver_templates}
  COMMAND perl -I ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/generateDrivers.pl
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Generating driver files"
  # VERBATIM
)
# If you decide you need to mess with DmrgDriver.pm you must reconfigure,
# building again is not sufficient

# for reasons that are not clear to me, almost every implementation file is
# referred to as a driver in the dmrg build, and many have non informative file
# names along that line don't make too much of it.

# I don't wrap all the libraries into a variable on purpose the build is so slow
# for such a small code and link time is substantial perhaps if we can reduce or
# chain the library dependencies better right now everything seems to need to
# link everything.

# Dmrg executable (main executable)
add_executable(dmrg $<TARGET_OBJECTS:Common> $<TARGET_OBJECTS:Su2Related> ${driver_templates} Qn.cpp dmrg.cpp)
target_compile_definitions(dmrg PUBLIC USE_BOOST)

target_link_libraries(dmrg PUBLIC kronutil psimaglite HDF5::HDF5 BLAS::BLAS LAPACK::LAPACK ${ADDITIONAL_LIBRARIES})
target_include_directories(dmrg PRIVATE ${CMAKE_CURRENT_BINARY_DIR} Engine)
target_sources(dmrg PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/GitRevision.h)

if(MPI_FOUND)
  target_compile_definitions(dmrg PUBLIC USE_MPI)
  target_link_libraries(dmrg PUBLIC MPI::MPI_CXX)
endif()

# DMRG++ as a library so cincuenta can use it
# We may instead bring cincuenta into this repo
add_library(dmrgpp STATIC $<TARGET_OBJECTS:Common> $<TARGET_OBJECTS:Su2Related> Qn.cpp ProgramGlobals.cpp)
target_include_directories(dmrgpp PRIVATE Engine)
target_link_libraries(dmrgpp PUBLIC psimaglite)

# Observe executable below
add_executable(
  observe
  observe.cpp
  ProgramGlobals.cpp
  Provenance.cpp
  Utils.cpp
  Su2Related.cpp
  Qn.cpp
  ObserveDriver0.cpp
  ObserveDriver1.cpp
  ObserveDriver2.cpp)

target_include_directories(observe PRIVATE ${CMAKE_CURRENT_BINARY_DIR} Engine)
target_compile_definitions(observe PUBLIC USE_BOOST)
target_sources(observe PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/GitRevision.h)
target_link_libraries(
  observe
  kronutil
  psimaglite::psimaglite
  HDF5::HDF5
  BLAS::BLAS
  LAPACK::LAPACK
  ${ADDITIONAL_LIBRARIES})

# Toolbox executable now
add_executable(toolboxdmrg toolboxdmrg.cpp Qn.cpp ProgramGlobals.cpp Provenance.cpp)
target_link_libraries(toolboxdmrg psimaglite::psimaglite HDF5::HDF5 ${ADDITIONAL_LIBRARIES})
target_sources(toolboxdmrg PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/GitRevision.h)
target_include_directories(toolboxdmrg PRIVATE ${CMAKE_CURRENT_BINARY_DIR} Engine)

# Many omegas executable
add_executable(omegas $<TARGET_OBJECTS:Common> $<TARGET_OBJECTS:Su2Related> Qn.cpp manyOmegas.cpp)
target_include_directories(omegas PRIVATE ${CMAKE_CURRENT_BINARY_DIR} Engine)
target_compile_definitions(omegas PUBLIC USE_BOOST)
target_link_libraries(
  omegas
  kronutil
  psimaglite::psimaglite
  HDF5::HDF5
  BLAS::BLAS
  LAPACK::LAPACK
  ${ADDITIONAL_LIBRARIES})

# Tests below for all executables
# Pipeline tests for all executables here

# For dmrg main executable
set(PATH_TO_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite/inputs)

add_test(NAME input2 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input2.ain)
# Enable this test after the bug fix applied to vsm
add_test(NAME input3 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input3.ain)
add_test(NAME input4 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input4.ain)
add_test(NAME input5 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input5.ain)
add_test(NAME input7 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input7.ain)
# input8 is too slow (disabled)
# add_test(NAME input8 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input8.ain)
add_test(NAME input10 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input10.ain)
# input11 and input12 need -DUSE_PTHREADS
# add_test(NAME input11 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input11.ain)
# add_test(NAME input12 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input12.ain)
add_test(NAME input20 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input20.ain)
add_test(NAME input21 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input21.ain)
add_test(NAME input30 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input30.ain)
add_test(NAME input32 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input32.ain)
add_test(NAME input70 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input70.ain)
add_test(NAME input71 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input71.ain)
add_test(NAME input120 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input120.ain)
add_test(NAME input121 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input121.ain  -p 12 "<gs|P0>")
add_test(NAME input122 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input122.ain  -p 12 "<gs|P0>,<gs|P1>,<P0|P1>")
add_test(NAME input123 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input123.ain  -p 12 "<gs|P0>,<P0|P1>,<gs|c'|P2>,<gs|c'|P3>,<P2|P3>")
add_test(NAME input126 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input126.ain)
add_test(NAME input127 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input127.ain  -p 12 "<gs|P0>")
add_test(NAME input170 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input170.ain)
# Enable this test after the bug fix applied to vsm
add_test(NAME input172 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input172.ain)
add_test(NAME input320 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input320.ain)
add_test(NAME input501 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input501.ain)
add_test(NAME input1300 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1300.ain)
add_test(NAME input1302 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1302.ain)
add_test(NAME input1304 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1304.ain)
add_test(NAME input1306 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1306.ain)
add_test(NAME input1814 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1814.ain)
add_test(NAME input1815 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1815.ain)
add_test(NAME input1817 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1817.ain)
# input1818 and 1819 --> too slow (disabled)
# add_test(NAME input1818 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1818.ain)
# add_test(NAME input1819 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input1819.ain)
add_test(NAME input2500 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input2500.ain)
add_test(NAME input2510 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input2510.ain)
# needs -DUSE_PTHREADS
# add_test(NAME input2700 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input2700.ain)
add_test(NAME input3020 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input3020.ain)
add_test(NAME input3021 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input3021.ain "<gs|c'|P2>")
add_test(NAME input3100 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input3100.ain)
# SDHS will remain disabled for a while
# add_test(NAME input4650 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input4650.ain)
add_test(NAME input4720 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input4720.ain)
add_test(NAME input6020 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input6020.ain)
add_test(NAME input6700 COMMAND ./dmrg -f ${PATH_TO_INPUTS}/input6700.ain)
