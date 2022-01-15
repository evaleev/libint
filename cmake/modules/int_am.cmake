# handle the defaulting and setting of the following variables
# * ENABLE_[ONEBODY|ERI2|ERI3|ERI|G12|G12DKH]
# * [LIBINT|ONEBODY|ERI2|ERI3|ERI|G12|G12DKH]_[MAX|OPT]_AM[|_LIST]
# * LIBINT_ONEBODY_DERIV
# * LIBINT_SUPPORTS_ONEBODY

# _candidate variables not needed for config.h but are used to figure
#   out the AM levels at the CMake level
#   so that libint2-config components may be defined and client codes can
#   require the detected library include gradient integrals up to AM=5 with
#   `find_package(Libint2 COMPONENTS twobody_c4_d1_l5)`


message(STATUS "LAB int.am")
message(STATUS "Processing integrals classes ...")

# <<<  overall derivatives level  >>>

set(_glob_classes_derivs ${ENABLE_ONEBODY};${ENABLE_ERI};${ENABLE_ERI3};${ENABLE_ERI2};${ENABLE_G12};${ENABLE_G12DKH})
list(SORT _glob_classes_derivs COMPARE NATURAL ORDER DESCENDING)  # CMake 3.18 for NATURAL
list(GET _glob_classes_derivs 0 _max_deriv)
message(STATUS "Preparing highest derivative level ${_max_deriv}")

# <<<  overall max_am defaults  >>>

list(LENGTH WITH_MAX_AM _ntokens_maxam)
if (_ntokens_maxam GREATER 1)
    math(EXPR _ntokens_xptd_max_deriv "${_max_deriv} + 1")
    if (NOT _ntokens_xptd_max_deriv EQUAL _ntokens_maxam)
        message(FATAL_ERROR "Invalid value for WITH_MAX_AM (${WITH_MAX_AM}). Highest ENABLE_ derivative (${_max_deriv}) requires list length ${_ntokens_xptd_max_deriv}, not ${_ntokens_maxam}.")
    endif()

    set(_sorted_WITH_MAX_AM ${WITH_MAX_AM})
    list(SORT _sorted_WITH_MAX_AM COMPARE NATURAL ORDER DESCENDING)
    list(GET _sorted_WITH_MAX_AM 0 _max_am)

    list(JOIN WITH_MAX_AM "," _sam)
    set(LIBINT_MAX_AM_LIST ${_sam})
    set(LIBINT_MAX_AM ${_max_am})  # only overall LIBINT, not specific integrals classes, sets both MAX_AM & MAX_AM_LIST
    set(_max_LIBINT_MAX_AM ${LIBINT_MAX_AM})
else()
    set(LIBINT_MAX_AM_LIST "")
    set(LIBINT_MAX_AM ${WITH_MAX_AM})
endif()

foreach(_d RANGE 0 ${_max_deriv})
    if (${_d} LESS _ntokens_maxam)
        list(GET WITH_MAX_AM ${_d} _candidate0_d${_d})
    else()
        set(_candidate0_d${_d} "-1")
    endif()
    message(VERBOSE "setting _candidate0_d${_d}=${_candidate0_d${_d}}")
endforeach()

if (LIBINT_MAX_AM GREATER_EQUAL 8)
    message(FATAL_ERROR "LIBINT_MAX_AM=${LIBINT_MAX_AM} is greater than 8. Are you sure you know what you are doing?")
elseif (LIBINT_MAX_AM LESS_EQUAL 0)
    message(FATAL_ERROR "Invalid value for LIBINT_MAX_AM (${LIBINT_MAX_AM}).")
endif()

message(STATUS "Preparing generic LIBINT_MAX_AM_LIST ${LIBINT_MAX_AM_LIST} and LIBINT_MAX_AM ${LIBINT_MAX_AM} for integrals class defaults.")

# <<<  overall opt_am defaults  >>>

list(LENGTH WITH_OPT_AM _ntokens_optam)
if (NOT WITH_OPT_AM EQUAL -1)
    if (NOT _ntokens_optam EQUAL _ntokens_maxam)
        # discard two cases: scalar opt and list max -and- list opt and scalar max
        message(FATAL_ERROR "Invalid format for WITH_OPT_AM (${WITH_OPT_AM}). Use the same format and length like `N` or `N0;N1;N2` as WITH_MAX_AM (${WITH_MAX_AM}).")
    endif()
endif()
if (_ntokens_optam GREATER 1)
    # list opt and list max: use list opt validating aginst max
    set(_processed_OPT_AM_LIST )
    math(EXPR _range_limit "${_ntokens_maxam} - 1")
    foreach(_d RANGE ${_range_limit})
        list(GET WITH_MAX_AM ${_d} _max_am)
        list(GET WITH_OPT_AM ${_d} _opt_am)
        if (_opt_am LESS_EQUAL _max_am)
            list(APPEND _processed_OPT_AM_LIST ${_opt_am})
        else()
            list(APPEND _processed_OPT_AM_LIST ${_max_am})
        endif()
    endforeach()

    list(JOIN _processed_OPT_AM_LIST "," LIBINT_OPT_AM_LIST)
    list(SORT _processed_OPT_AM_LIST COMPARE NATURAL ORDER DESCENDING)
    list(GET _processed_OPT_AM_LIST 0 LIBINT_OPT_AM)
else()
    if(WITH_OPT_AM EQUAL -1)
        if (_ntokens_maxam GREATER 1)
            # no opt and list max: default list opt from max
            set(_processed_OPT_AM_LIST )
            math(EXPR _range_limit "${_ntokens_maxam} - 1")
            foreach(_d RANGE ${_range_limit})
                list(GET WITH_MAX_AM ${_d} _max_am)
                math(EXPR _opt_am "${_max_am}/2 + 1")
                list(APPEND _processed_OPT_AM_LIST ${_opt_am})
            endforeach()

            list(JOIN _processed_OPT_AM_LIST "," LIBINT_OPT_AM_LIST)
            list(SORT _processed_OPT_AM_LIST COMPARE NATURAL ORDER DESCENDING)
            list(GET _processed_OPT_AM_LIST 0 LIBINT_OPT_AM)
        else()
            # no opt and scalar max: default scalar opt from max
            set(LIBINT_OPT_AM_LIST "")
            math(EXPR LIBINT_OPT_AM "${LIBINT_MAX_AM}/2 + 1")
        endif()
    else()
        # scalar opt and scalar max: use scalar opt validating aginst max
        set(LIBINT_OPT_AM_LIST "")
        set(LIBINT_OPT_AM ${WITH_OPT_AM})

        if (LIBINT_OPT_AM GREATER LIBINT_MAX_AM)
            set(LIBINT_OPT_AM ${LIBINT_MAX_AM})
        endif()
    endif()
endif()

message(STATUS "Preparing generic LIBINT_OPT_AM_LIST ${LIBINT_OPT_AM_LIST} and LIBINT_OPT_AM ${LIBINT_OPT_AM} for integrals class defaults.")

# <<<  Macro  >>>

macro(process_integrals_class class)

    list(LENGTH ENABLE_${class} _ntokens)
    if (NOT _ntokens EQUAL 1)
        message(FATAL_ERROR "Invalid value for ENABLE_${class} (${ENABLE_${class}}). Use scalar of maximum derivative level, not list.")
    endif()

    if (ENABLE_${class} GREATER_EQUAL 0)
        set(INCLUDE_${class} ${ENABLE_${class}})

        foreach(_d RANGE 0 ${_max_deriv})
            if (${_d} LESS_EQUAL ${INCLUDE_${class}})
                set(_candidate0_${class}_d${_d} ${_candidate0_d${_d}})
                message(VERBOSE "setting _candidate0_${class}_d${_d}=${_candidate0_${class}_d${_d}}")
            endif()
        endforeach()

        set(LIBINT_SUPPORTS_${class} yes)
        set(LIBINT_${class}_DERIV ${INCLUDE_${class}})
        message(STATUS "Enabling integrals class ${class} to derivative ${INCLUDE_${class}}")
    else()
        set(INCLUDE_${class} "-1")
        set(${class}_MAX_AM "")
        set(${class}_MAX_AM_LIST "")
        message(STATUS "Disabling integrals class ${class}")
    endif()

    if (ENABLE_${class} GREATER_EQUAL 0)
        list(LENGTH WITH_${class}_MAX_AM _ntokens)
        if (_ntokens GREATER 1)
            math(EXPR _ntokens_xptd_max_deriv "${INCLUDE_${class}} + 1")
            if (NOT _ntokens_xptd_max_deriv EQUAL _ntokens)
                message(FATAL_ERROR "Invalid value for WITH_${class}_MAX_AM (${WITH_${class}_MAX_AM}). ENABLE_${class} derivative (${INCLUDE_${class}}) requires list length ${_ntokens_xptd_max_deriv}, not ${_ntokens}.")
            endif()

            foreach(_d RANGE ${INCLUDE_${class}})
                list(GET WITH_${class}_MAX_AM ${_d} _candidate_${class}_d${_d})
                message(VERBOSE "setting _candidate_${class}_d${_d}=${_candidate_${class}_d${_d}}")

                if (_candidate_${class}_d${_d} GREATER LIBINT_MAX_AM)
                    message(FATAL_ERROR "Invalid value for WITH_${class}_MAX_AM derivative element ${_d} (${_candidate_${class}_d${_d}} > ${LIBINT_MAX_AM}).")
                    # note this check is necessary but insufficient since per-d max may be available from LIBINT_MAX_AM_LIST.
                    # suggest requiring all WITH_*_AM options to be either `N` or `N0;N1;N2` format for cleaner validation.
                elseif (_candidate_${class}_d${_d} LESS_EQUAL 0)
                    message(FATAL_ERROR "Invalid value for WITH_${class}_MAX_AM derivative element ${_d} (${_candidate_${class}_d${_d}} <= 0).")
                endif()
            endforeach()

            list(JOIN WITH_${class}_MAX_AM "," ${class}_MAX_AM_LIST)
            set(${class}_MAX_AM "")
        else()
            set(${class}_MAX_AM_LIST "")
            if (WITH_${class}_MAX_AM EQUAL -1)
                foreach(_d RANGE ${INCLUDE_${class}})
                    if (${_candidate0_${class}_d${_d}} EQUAL -1)
                        set(_candidate_${class}_d${_d} ${_candidate0_${class}_d0})
                    else()
                        set(_candidate_${class}_d${_d} ${_candidate0_${class}_d${_d}})
                    endif()
                    message(VERBOSE "setting _candidate_${class}_d${_d}=${_candidate_${class}_d${_d}}")
                endforeach()

                set(${class}_MAX_AM "")
                # note: could set class_MAX_AM/LIST from default (in configure.ac, looks like at least scalar var set)
                #       but philosophy is to set user-only intent and leave further defaulting to compiled code. wrong?
            else()
                set(${class}_MAX_AM ${WITH_${class}_MAX_AM})

                foreach(_d RANGE ${INCLUDE_${class}})
                    set(_candidate_${class}_d${_d} ${${class}_MAX_AM})
                    message(VERBOSE "setting _candidate_${class}_d${_d}=${_candidate_${class}_d${_d}}")
                endforeach()

                if (${class}_MAX_AM GREATER_EQUAL 8)
                    message(FATAL_ERROR "Value for ${class}_MAX_AM too high (${${class}_MAX_AM} >= 8). Are you sure you know what you are doing?")
                elseif (${class}_MAX_AM LESS_EQUAL 0)
                    message(FATAL_ERROR "Invalid value for ${class}_MAX_AM (${${class}_MAX_AM} <= 0).")
                endif()
            endif()
        endif()
        if (LIBINT_MAX_AM_LIST)
            set(_msg ${LIBINT_MAX_AM_LIST})
        else()
            set(_msg ${LIBINT_MAX_AM})
        endif()
        message(STATUS "Enabling integrals class ${class} to max AM ${${class}_MAX_AM}${${class}_MAX_AM_LIST} (else ${_msg})")

        list(LENGTH WITH_${class}_OPT_AM _ntokens)
        if (_ntokens GREATER 1)
            if (NOT _ntokens_xptd_max_deriv EQUAL _ntokens)
                message(FATAL_ERROR "Invalid value for WITH_${class}_OPT_AM (${WITH_${class}_OPT_AM}). ENABLE_${class} derivative (${INCLUDE_${class}}) requires list length ${_ntokens_xptd_max_deriv}, not ${_ntokens}.")
            endif()

            list(JOIN WITH_${class}_OPT_AM "," ${class}_OPT_AM_LIST)
            set(${class}_OPT_AM "")
        else()
            set(${class}_OPT_AM_LIST "")
            if (WITH_${class}_OPT_AM EQUAL -1)
                set(${class}_OPT_AM "")
            else()
                set(${class}_OPT_AM ${WITH_${class}_OPT_AM})
            endif()
        endif()
        if (LIBINT_OPT_AM_LIST)
            set(_msg ${LIBINT_OPT_AM_LIST})
        else()
            set(_msg ${LIBINT_OPT_AM})
        endif()
        message(STATUS "Enabling integrals class ${class} to opt AM ${${class}_OPT_AM}${${class}_OPT_AM_LIST} (else ${_msg})")
    endif()
endmacro()


process_integrals_class(ONEBODY)
process_integrals_class(ERI)
process_integrals_class(ERI3)
process_integrals_class(ERI2)

# discrepancy, as configure doesn't do AM_LIST for these
process_integrals_class(G12)
process_integrals_class(G12DKH)

# form list of active <class>_<deriv><max_am> strings to use in Libint2Config
set(Libint2_ERI_COMPONENTS "")
foreach(_cls ERI;ERI3;ERI2)
    if (_cls STREQUAL "ERI")
        set(_lbl "eri_c4")
    elseif (_cls STREQUAL "ERI3")
        set(_lbl "eri_c3")
    elseif (_cls STREQUAL "ERI2")
        set(_lbl "eri_c2")
    endif()

    if (INCLUDE_${_cls} GREATER -1)
        foreach (_d RANGE 0 ${INCLUDE_${_cls}})
            foreach(_l RANGE 2 ${_candidate_${_cls}_d${_d}})
                list(APPEND Libint2_ERI_COMPONENTS "${_lbl}_d${_d}_l${_l}")
                message(VERBOSE "setting component ${_lbl}_d${_d}_l${_l}")
            endforeach()
        endforeach()
    endif()
endforeach()
message(STATUS "Library will satisfy ERI AM components: ${Libint2_ERI_COMPONENTS}")

