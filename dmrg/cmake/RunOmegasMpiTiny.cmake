if(NOT DEFINED DMRG_EXECUTABLE OR NOT EXISTS "${DMRG_EXECUTABLE}")
  message(FATAL_ERROR "DMRG_EXECUTABLE is not set or does not exist: ${DMRG_EXECUTABLE}")
endif()

if(NOT DEFINED OMEGAS_EXECUTABLE OR NOT EXISTS "${OMEGAS_EXECUTABLE}")
  message(FATAL_ERROR "OMEGAS_EXECUTABLE is not set or does not exist: ${OMEGAS_EXECUTABLE}")
endif()

if(NOT DEFINED GS_INPUT OR NOT EXISTS "${GS_INPUT}")
  message(FATAL_ERROR "GS_INPUT is not set or does not exist: ${GS_INPUT}")
endif()

if(NOT DEFINED OMEGAS_INPUT OR NOT EXISTS "${OMEGAS_INPUT}")
  message(FATAL_ERROR "OMEGAS_INPUT is not set or does not exist: ${OMEGAS_INPUT}")
endif()

if(NOT DEFINED TEST_WORKDIR OR TEST_WORKDIR STREQUAL "")
  message(FATAL_ERROR "TEST_WORKDIR is not set")
endif()

if(NOT DEFINED MPIEXEC_EXECUTABLE OR MPIEXEC_EXECUTABLE STREQUAL "")
  message(FATAL_ERROR "MPIEXEC_EXECUTABLE is not set")
endif()

if(NOT DEFINED MPIEXEC_NUMPROC_FLAG OR MPIEXEC_NUMPROC_FLAG STREQUAL "")
  set(MPIEXEC_NUMPROC_FLAG "-n")
endif()

file(REMOVE_RECURSE "${TEST_WORKDIR}")
file(MAKE_DIRECTORY "${TEST_WORKDIR}")

execute_process(
  COMMAND "${DMRG_EXECUTABLE}" -f "${GS_INPUT}"
  WORKING_DIRECTORY "${TEST_WORKDIR}"
  RESULT_VARIABLE _gs_result
  OUTPUT_VARIABLE _gs_stdout
  ERROR_VARIABLE _gs_stderr)

set(_gs_output "${_gs_stdout}\n${_gs_stderr}")
if(NOT _gs_result EQUAL 0)
  message(FATAL_ERROR "tiny GS run failed with exit code ${_gs_result}\n${_gs_output}")
endif()

if(NOT EXISTS "${TEST_WORKDIR}/tinyGs.hd5")
  message(FATAL_ERROR "tiny GS run did not create ${TEST_WORKDIR}/tinyGs.hd5\n${_gs_output}")
endif()

set(_mpi_cmd
    "${MPIEXEC_EXECUTABLE}"
    "${MPIEXEC_NUMPROC_FLAG}" "2")

if(DEFINED MPIEXEC_PREFLAGS AND NOT MPIEXEC_PREFLAGS STREQUAL "")
  list(APPEND _mpi_cmd ${MPIEXEC_PREFLAGS})
endif()

list(APPEND _mpi_cmd "${OMEGAS_EXECUTABLE}")

if(DEFINED MPIEXEC_POSTFLAGS AND NOT MPIEXEC_POSTFLAGS STREQUAL "")
  list(APPEND _mpi_cmd ${MPIEXEC_POSTFLAGS})
endif()

list(APPEND _mpi_cmd -f "${OMEGAS_INPUT}" -O tinyOmega)

execute_process(
  COMMAND ${_mpi_cmd}
  WORKING_DIRECTORY "${TEST_WORKDIR}"
  RESULT_VARIABLE _omegas_result
  OUTPUT_VARIABLE _omegas_stdout
  ERROR_VARIABLE _omegas_stderr)

set(_omegas_output "${_omegas_stdout}\n${_omegas_stderr}")
if(NOT _omegas_result EQUAL 0)
  message(FATAL_ERROR "tiny omegas MPI run failed with exit code ${_omegas_result}\n${_omegas_output}")
endif()

foreach(_rank IN ITEMS 0 1)
  string(FIND "${_omegas_output}" "MPI rank=${_rank}" _rank_seen)
  if(_rank_seen EQUAL -1)
    message(FATAL_ERROR "tiny omegas MPI run did not show work on rank ${_rank}\n${_omegas_output}")
  endif()
endforeach()

foreach(_idx IN ITEMS 0 1)
  if(NOT EXISTS "${TEST_WORKDIR}/tinyOmega${_idx}.hd5")
    message(FATAL_ERROR "tiny omegas MPI run did not create ${TEST_WORKDIR}/tinyOmega${_idx}.hd5\n${_omegas_output}")
  endif()
endforeach()

message(STATUS "tiny omegas MPI run completed and produced both omega outputs")
