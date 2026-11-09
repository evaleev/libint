# Checks that the generated code does not depend on the order in which the
# shell sets and the set-level recurrence relations (RRs) are generated (which
# makes it independent of the number of processes that generate it, see
# LIBINT2_GENERATOR_JOBS): runs the generator serially with the shell sets and
# the RRs in other orders, and compares its output with the regular output.
#
# Usage: cmake -DBUILD_LIBINT=<generator> -DREFERENCE_DIR=<regular output dir>
#              -DWORK_DIR=<scratch dir> [-DORDERS=<order>[,<order>...]]
#              -P check_codegen_order.cmake
# where <order> is <rr order>/<shell set order>, the values of LIBINT_RR_ORDER
# and LIBINT_SHELLSET_ORDER (reverse, label, shuffle[:<seed>]). Different
# orders expose different defects, hence the default checks two: in a broad
# configuration (1-body derivatives with multipoles, ERI, ERI3, ERI2, RKB,
# G12; MAX_AM=2) the RR order label exposed a precomputed RR argument leaking
# into the code and shuffle:5 the numbering of automatic variables in stack
# order, which reverse and shuffle:1..4 did not expose. The shuffle is the
# same on every platform (see permute_for_testing() in
# src/bin/libint/codegen_jobs.h), hence so is the order it checks.

foreach(_var BUILD_LIBINT REFERENCE_DIR WORK_DIR)
  if (NOT DEFINED ${_var})
    message(FATAL_ERROR "check_codegen_order.cmake: ${_var} must be defined")
  endif()
endforeach()
if (NOT DEFINED ORDERS)
  set(ORDERS "label/reverse,shuffle:5/shuffle:1")
endif()
string(REPLACE "," ";" ORDERS "${ORDERS}")

set(_failed_orders)
foreach(ORDER IN LISTS ORDERS)
  string(REPLACE "/" ";" _orders "${ORDER}")
  list(GET _orders 0 _rr_order)
  list(GET _orders 1 _shellset_order)
  # '/' and ':' are not valid in (Windows) file names
  string(REGEX REPLACE "[/:]" "-" _order_tag "${ORDER}")
  file(REMOVE_RECURSE "${WORK_DIR}")
  file(MAKE_DIRECTORY "${WORK_DIR}")
  execute_process(
    # serially, in this order; the reference may have been generated in
    # parallel (LIBINT2_GENERATOR_JOBS), hence this also checks that parallel
    # generation does not change the generated code
    COMMAND ${CMAKE_COMMAND} -E env LIBINT_RR_ORDER=${_rr_order}
            LIBINT_SHELLSET_ORDER=${_shellset_order} LIBINT_GENERATOR_JOBS=1
            "${BUILD_LIBINT}"
    WORKING_DIRECTORY "${WORK_DIR}"
    OUTPUT_FILE "${WORK_DIR}-${_order_tag}.log"
    ERROR_FILE "${WORK_DIR}-${_order_tag}.log"
    RESULT_VARIABLE _result)
  if (NOT _result EQUAL 0)
    message(FATAL_ERROR "generator failed with order ${ORDER} (${_result}), see ${WORK_DIR}-${_order_tag}.log")
  endif()

  # every file must be identical to the regular output
  file(GLOB _files RELATIVE "${WORK_DIR}" "${WORK_DIR}/*")
  list(LENGTH _files _nfiles)
  set(_ndiffer 0)
  # ... which must not contain more files (e.g. left behind by the workers
  # that generated it), except stale ones from the runs of earlier generators
  # (older than this one)
  file(GLOB _reference_files LIST_DIRECTORIES true RELATIVE "${REFERENCE_DIR}"
       "${REFERENCE_DIR}/*" "${REFERENCE_DIR}/.libint*")
  list(REMOVE_DUPLICATES _reference_files)
  foreach(_f IN LISTS _reference_files)
    if (NOT EXISTS "${WORK_DIR}/${_f}" AND
        "${REFERENCE_DIR}/${_f}" IS_NEWER_THAN "${BUILD_LIBINT}")
      message(SEND_ERROR "${_f}: generated in the regular output only")
      math(EXPR _ndiffer "${_ndiffer} + 1")
    endif()
  endforeach()
  foreach(_f IN LISTS _files)
    if (NOT EXISTS "${REFERENCE_DIR}/${_f}")
      message(SEND_ERROR "${_f}: generated with order ${ORDER} only")
      math(EXPR _ndiffer "${_ndiffer} + 1")
      continue()
    endif()
    file(SHA256 "${WORK_DIR}/${_f}" _h1)
    file(SHA256 "${REFERENCE_DIR}/${_f}" _h2)
    if (NOT _h1 STREQUAL _h2)
      message(SEND_ERROR "${_f}: differs with order ${ORDER}")
      math(EXPR _ndiffer "${_ndiffer} + 1")
    endif()
  endforeach()
  if (_ndiffer GREATER 0)
    message(SEND_ERROR "${_ndiffer} of ${_nfiles} generated files depend on the order of the shell sets and the RRs (order ${ORDER}); the output is kept in ${WORK_DIR}-${_order_tag}")
    list(APPEND _failed_orders "${ORDER}")
    file(REMOVE_RECURSE "${WORK_DIR}-${_order_tag}")
    file(RENAME "${WORK_DIR}" "${WORK_DIR}-${_order_tag}")
  else()
    message(STATUS "all ${_nfiles} generated files are independent of the order of the shell sets and the RRs (order ${ORDER})")
    file(REMOVE_RECURSE "${WORK_DIR}")
  endif()
endforeach()

if (_failed_orders)
  message(FATAL_ERROR "the generated code depends on the order of the shell sets and the RRs (failed: ${_failed_orders})")
endif()
