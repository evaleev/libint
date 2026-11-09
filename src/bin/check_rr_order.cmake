# Checks that the generated code does not depend on the order in which the
# set-level recurrence relations (RRs) are generated: runs the generator with
# the RRs in other orders and compares its output with the regular output.
#
# Usage: cmake -DBUILD_LIBINT=<generator> -DREFERENCE_DIR=<regular output dir>
#              -DWORK_DIR=<scratch dir> [-DORDERS=<order>[,<order>...]]
#              -P check_rr_order.cmake
# where <order> is a value of LIBINT_RR_ORDER (reverse, label, shuffle[:<seed>]).
# Each order exposes different defects (e.g. label and shuffle:2 each exposed
# one that the other orders did not), hence the default checks two.

foreach(_var BUILD_LIBINT REFERENCE_DIR WORK_DIR)
  if (NOT DEFINED ${_var})
    message(FATAL_ERROR "check_rr_order.cmake: ${_var} must be defined")
  endif()
endforeach()
if (NOT DEFINED ORDERS)
  set(ORDERS "label,shuffle:2")
endif()
string(REPLACE "," ";" ORDERS "${ORDERS}")

set(_failed_orders)
foreach(ORDER IN LISTS ORDERS)
  file(REMOVE_RECURSE "${WORK_DIR}")
  file(MAKE_DIRECTORY "${WORK_DIR}")
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E env LIBINT_RR_ORDER=${ORDER} "${BUILD_LIBINT}"
    WORKING_DIRECTORY "${WORK_DIR}"
    OUTPUT_FILE "${WORK_DIR}-${ORDER}.log"
    ERROR_FILE "${WORK_DIR}-${ORDER}.log"
    RESULT_VARIABLE _result)
  if (NOT _result EQUAL 0)
    message(FATAL_ERROR "generator failed with LIBINT_RR_ORDER=${ORDER} (${_result}), see ${WORK_DIR}-${ORDER}.log")
  endif()

  # every file must be identical to the regular output (which may contain stale
  # files from earlier generator runs, hence compare in this direction only)
  file(GLOB _files RELATIVE "${WORK_DIR}" "${WORK_DIR}/*")
  list(LENGTH _files _nfiles)
  set(_ndiffer 0)
  foreach(_f IN LISTS _files)
    if (NOT EXISTS "${REFERENCE_DIR}/${_f}")
      message(SEND_ERROR "${_f}: generated with LIBINT_RR_ORDER=${ORDER} only")
      math(EXPR _ndiffer "${_ndiffer} + 1")
      continue()
    endif()
    file(SHA256 "${WORK_DIR}/${_f}" _h1)
    file(SHA256 "${REFERENCE_DIR}/${_f}" _h2)
    if (NOT _h1 STREQUAL _h2)
      message(SEND_ERROR "${_f}: differs with LIBINT_RR_ORDER=${ORDER}")
      math(EXPR _ndiffer "${_ndiffer} + 1")
    endif()
  endforeach()
  if (_ndiffer GREATER 0)
    message(SEND_ERROR "${_ndiffer} of ${_nfiles} generated files depend on the order of the RRs (LIBINT_RR_ORDER=${ORDER}); the output is kept in ${WORK_DIR}-${ORDER}")
    list(APPEND _failed_orders "${ORDER}")
    file(REMOVE_RECURSE "${WORK_DIR}-${ORDER}")
    file(RENAME "${WORK_DIR}" "${WORK_DIR}-${ORDER}")
  else()
    message(STATUS "all ${_nfiles} generated files are independent of the order of the RRs (LIBINT_RR_ORDER=${ORDER})")
    file(REMOVE_RECURSE "${WORK_DIR}")
  endif()
endforeach()

if (_failed_orders)
  message(FATAL_ERROR "the generated code depends on the order of the RRs (failed: ${_failed_orders})")
endif()
