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
