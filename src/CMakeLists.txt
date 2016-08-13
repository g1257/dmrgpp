cmake_minimum_required(VERSION 3.4)
project (dmrgpp)

cmake_policy(SET CMP0056 NEW)
cmake_policy(SET CMP0065 NEW)
cmake_policy(SET CMP0042 NEW)

include_directories(${PsimagLite_SOURCE_DIR}
  ${PsimagLite_SOURCE_DIR}/../
  Engine
  )

add_library( Common OBJECT ProgramGlobals.cpp 
  Provenance.cpp 
  Utils.cpp
  )
add_library( Su2Related OBJECT Su2Related.cpp )

set(driver_list DmrgDriver0.cpp DmrgDriver1.cpp DmrgDriver2.cpp DmrgDriver3.cpp DmrgDriver4.cpp DmrgDriver5.cpp DmrgDriver6.cpp DmrgDriver7.cpp DmrgDriver8.cpp DmrgDriver9.cpp DmrgDriver10.cpp DmrgDriver11.cpp DmrgDriver12.cpp DmrgDriver13.cpp DmrgDriver14.cpp DmrgDriver15.cpp DmrgDriver16.cpp DmrgDriver17.cpp DmrgDriver18.cpp DmrgDriver19.cpp DmrgDriver20.cpp DmrgDriver21.cpp DmrgDriver22.cpp DmrgDriver23.cpp)

foreach(driver ${driver_list})
  list(APPEND driver_templates ${CMAKE_CURRENT_SOURCE_DIR}/${driver})
endforeach()

message("${driver_templates}")
add_custom_command( OUTPUT ${driver_templates}
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_drivers.pl
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Generating driver files"
#  VERBATIM
  )

add_executable (dmrg $<TARGET_OBJECTS:Common> $<TARGET_OBJECTS:Su2Related>
  RestartStruct.cpp FiniteLoop.cpp DmrgDriver0.cpp DmrgDriver1.cpp DmrgDriver2.cpp DmrgDriver3.cpp DmrgDriver4.cpp DmrgDriver5.cpp DmrgDriver6.cpp DmrgDriver7.cpp DmrgDriver8.cpp DmrgDriver9.cpp DmrgDriver10.cpp DmrgDriver11.cpp DmrgDriver12.cpp DmrgDriver13.cpp DmrgDriver14.cpp DmrgDriver15.cpp DmrgDriver16.cpp DmrgDriver17.cpp DmrgDriver18.cpp DmrgDriver19.cpp DmrgDriver20.cpp DmrgDriver21.cpp DmrgDriver22.cpp DmrgDriver23.cpp dmrg.cpp)

add_executable (toolboxdmrg $<TARGET_OBJECTS:Common> toolboxdmrg.cpp)

add_executable (observe $<TARGET_OBJECTS:Common> $<TARGET_OBJECTS:Su2Related> observe.cpp)

target_link_libraries (dmrg
  PsimagLite
  ${ADDITIONAL_LIBRARIES}
  )

target_link_libraries (toolboxdmrg
  PsimagLite
  ${ADDITIONAL_LIBRARIES}
  )

target_link_libraries (observe
  PsimagLite
  ${ADDITIONAL_LIBRARIES}
  )

