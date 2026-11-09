# Checks how the generator handles a worker process that dies (e.g. killed by
# the out-of-memory killer) with LIBINT_GENERATOR_JOBS > 1: the generator must
# fail with an exit code (not die by a signal, e.g. SIGPIPE), say how the worker
# died, and remove the directories of the workers, which hold partial output
# (and would otherwise be exported with the library, see populate.cmake).
#
# Usage: cmake -DBUILD_LIBINT=<generator> -DWORK_DIR=<scratch dir>
#              -P check_codegen_failure.cmake
#
# LIBINT_GENERATOR_TEST_KILL_WORKER=<k>:<point> (see codegen_jobs.cc) kills
# worker k with SIGKILL at <point>:
#   shellsets: while the workers generate the shell sets
#   reported: after the worker reported its shell sets, i.e. when the generator
#             sends it the RRs to generate (a write to a closed pipe)

foreach(_var BUILD_LIBINT WORK_DIR)
  if (NOT DEFINED ${_var})
    message(FATAL_ERROR "check_codegen_failure.cmake: ${_var} must be defined")
  endif()
endforeach()

set(_failed)
foreach(_point shellsets reported)
  file(REMOVE_RECURSE "${WORK_DIR}")
  file(MAKE_DIRECTORY "${WORK_DIR}")
  set(_log "${WORK_DIR}-${_point}.log")
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E env LIBINT_GENERATOR_JOBS=2
            LIBINT_GENERATOR_TEST_KILL_WORKER=1:${_point}
            "${BUILD_LIBINT}"
    WORKING_DIRECTORY "${WORK_DIR}"
    OUTPUT_FILE "${_log}"
    ERROR_FILE "${_log}"
    RESULT_VARIABLE _result)
  file(READ "${_log}" _output)

  if (NOT _result MATCHES "^[1-9][0-9]*$")
    message(SEND_ERROR "${_point}: the generator should fail with an exit code, got \"${_result}\"; see ${_log}")
    list(APPEND _failed ${_point})
  endif()
  if (NOT _output MATCHES "process 1 generating the [a-zA-Z ]+ was killed by signal")
    message(SEND_ERROR "${_point}: the generator did not report how the worker died; see ${_log}")
    list(APPEND _failed ${_point})
  endif()
  file(GLOB _left RELATIVE "${WORK_DIR}" "${WORK_DIR}/.libint_worker*" "${WORK_DIR}/*.libint-rename")
  if (_left)
    message(SEND_ERROR "${_point}: the generator left ${_left} behind in ${WORK_DIR}")
    list(APPEND _failed ${_point})
  endif()

  if (NOT _point IN_LIST _failed)
    message(STATUS "${_point}: the generator failed cleanly")
    file(REMOVE "${_log}")
  endif()
endforeach()
file(REMOVE_RECURSE "${WORK_DIR}")

if (_failed)
  message(FATAL_ERROR "the generator mishandles failed workers (${_failed})")
endif()
