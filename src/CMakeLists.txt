cmake_minimum_required(VERSION 3.21)

find_package(PsimagLite REQUIRED CONFIG)

add_custom_command(
  OUTPUT GitRevision.h
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/createGitRevision.pl GitRevision.h
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/createGitRevision.pl
  COMMENT "Generating Git Revision File"
  VERBATIM)

add_subdirectory(KronUtil)

add_library(Common OBJECT ProgramGlobals.cpp Provenance.cpp Utils.cpp)
target_sources(Common PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/GitRevision.h)
target_link_libraries(Common PRIVATE psimaglite::psimaglite psimaglite::loki
                                     HDF5::HDF5)
target_include_directories(Common PRIVATE ${CMAKE_CURRENT_BINARY_DIR} Engine)

add_library(Su2Related OBJECT Su2Related.cpp)
target_link_libraries(Su2Related PRIVATE psimaglite::psimaglite HDF5::HDF5)
target_include_directories(Su2Related PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
                                              Engine)

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
  COMMAND perl -I ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/generateDrivers.pl
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Generating driver files"
  # VERBATIM
)
# If you decide you need to mess with DmrgDriver.pm you must reconfigure,
# building again is not sufficient

# for reasons that are not clear to me, almost every implementation file is
# referred to as a driver in the dmrg build, and many have non informative file
# names along that line don't make too much of it.
add_library(toolboxdmrg $<TARGET_OBJECTS:Common> toolboxdmrg.cpp)
target_link_libraries(toolboxdmrg psimaglite::psimaglite HDF5::HDF5
                      ${ADDITIONAL_LIBRARIES})
target_sources(toolboxdmrg PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/GitRevision.h)
target_include_directories(toolboxdmrg PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
                                               Engine)

set(ObserveDriver0.cpp ObserveDriver1.cpp ObserveDriver2.cpp)
add_library(observe $<TARGET_OBJECTS:Common> $<TARGET_OBJECTS:Su2Related>
                    ${observe_drivers} observe.cpp)
target_include_directories(observe PRIVATE ${CMAKE_CURRENT_BINARY_DIR} Engine)
target_link_libraries(observe psimaglite::psimaglite HDF5::HDF5 BLAS::BLAS
                      ${ADDITIONAL_LIBRARIES})

add_library(omegas manyOmegas.cpp procOmegas.cpp)
target_include_directories(omegas PRIVATE ${CMAKE_CURRENT_BINARY_DIR} Engine)
target_link_libraries(omegas psimaglite::psimaglite HDF5::HDF5)

# I don't wrap all the libraries into a variable on purpose the build is so slow
# for such a small code and link time is substantial perhaps if we can reduce or
# chain the library dependencies better right now everything seems to need to
# link everything.
add_executable(dmrg $<TARGET_OBJECTS:Common> $<TARGET_OBJECTS:Su2Related>
                    ${driver_templates} Qn.cpp dmrg.cpp)
target_link_libraries(
  dmrg
  toolboxdmrg
  observe
  omegas
  kronutil
  psimaglite::psimaglite
  HDF5::HDF5
  BLAS::BLAS
  ${ADDITIONAL_LIBRARIES})
target_include_directories(dmrg PRIVATE ${CMAKE_CURRENT_BINARY_DIR} Engine)
target_sources(dmrg PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/GitRevision.h)
