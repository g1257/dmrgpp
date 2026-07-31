if(NOT DEFINED OMEGAS_EXECUTABLE OR NOT EXISTS "${OMEGAS_EXECUTABLE}")
  message(FATAL_ERROR "OMEGAS_EXECUTABLE is not set or does not exist: ${OMEGAS_EXECUTABLE}")
endif()

if(NOT DEFINED OMEGAS_INPUT OR NOT EXISTS "${OMEGAS_INPUT}")
  message(FATAL_ERROR "OMEGAS_INPUT is not set or does not exist: ${OMEGAS_INPUT}")
endif()

if(NOT DEFINED MPIEXEC_EXECUTABLE OR MPIEXEC_EXECUTABLE STREQUAL "")
  message(FATAL_ERROR "MPIEXEC_EXECUTABLE is not set")
endif()

if(NOT DEFINED MPIEXEC_NUMPROC_FLAG OR MPIEXEC_NUMPROC_FLAG STREQUAL "")
  set(MPIEXEC_NUMPROC_FLAG "-n")
endif()

set(_cmd
    "${MPIEXEC_EXECUTABLE}"
    "${MPIEXEC_NUMPROC_FLAG}" "2")

if(DEFINED MPIEXEC_PREFLAGS AND NOT MPIEXEC_PREFLAGS STREQUAL "")
  list(APPEND _cmd ${MPIEXEC_PREFLAGS})
endif()

list(APPEND _cmd "${OMEGAS_EXECUTABLE}")

if(DEFINED MPIEXEC_POSTFLAGS AND NOT MPIEXEC_POSTFLAGS STREQUAL "")
  list(APPEND _cmd ${MPIEXEC_POSTFLAGS})
endif()

list(APPEND _cmd -d -f "${OMEGAS_INPUT}")

execute_process(
  COMMAND ${_cmd}
  RESULT_VARIABLE _result
  OUTPUT_VARIABLE _stdout
  ERROR_VARIABLE _stderr)

set(_output "${_stdout}\n${_stderr}")

if(NOT _result EQUAL 0)
  message(FATAL_ERROR "omegas MPI dry run failed with exit code ${_result}\n${_output}")
endif()

string(FIND "${_output}" "MPI rank=0" _rank0)
string(FIND "${_output}" "MPI rank=1" _rank1)
string(FIND "${_output}" "NOT done because -d" _dryrun)

if(_rank0 EQUAL -1 OR _rank1 EQUAL -1)
  message(FATAL_ERROR "omegas MPI dry run did not show work on both ranks\n${_output}")
endif()

if(_dryrun EQUAL -1)
  message(FATAL_ERROR "omegas MPI dry run did not exercise dry-run path\n${_output}")
endif()

message(STATUS "omegas MPI dry run exercised both ranks")
