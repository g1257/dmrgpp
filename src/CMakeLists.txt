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

set(driver_list
    DmrgDriver0.cpp
    DmrgDriver1.cpp
    DmrgDriver2.cpp
    DmrgDriver3.cpp
    DmrgDriver4.cpp
    DmrgDriver5.cpp
    DmrgDriver6.cpp
    DmrgDriver7.cpp
    DmrgDriver8.cpp
    DmrgDriver9.cpp
    DmrgDriver10.cpp
    DmrgDriver11.cpp
    DmrgDriver12.cpp
    DmrgDriver13.cpp
    DmrgDriver14.cpp
    DmrgDriver15.cpp
    DmrgDriver16.cpp
    DmrgDriver17.cpp
    DmrgDriver18.cpp
    DmrgDriver19.cpp
    DmrgDriver20.cpp
    DmrgDriver21.cpp
    DmrgDriver22.cpp
    DmrgDriver23.cpp)

foreach(driver ${driver_list})
  list(APPEND driver_templates ${CMAKE_CURRENT_SOURCE_DIR}/${driver})
endforeach()

message("${driver_templates}")
add_custom_command(
  OUTPUT ${driver_templates}
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/generateDrivers.pl
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Generating driver files"
  # VERBATIM
)

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

add_executable(dmrg $<TARGET_OBJECTS:Common> $<TARGET_OBJECTS:Su2Related>
                    ${driver_list} Qn.cpp dmrg.cpp)
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
