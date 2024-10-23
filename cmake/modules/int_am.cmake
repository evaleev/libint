# handle the defaulting and setting of the following variables
# * ENABLE_[ONEBODY|ERI2|ERI3|ERI|G12|G12DKH]
# * [LIBINT|ONEBODY|ERI2|ERI3|ERI|G12|G12DKH]_[MAX|OPT]_AM[|_LIST]
# * MULTIPOLE_MAX_ORDER
# * LIBINT_ONEBODY_DERIV
# * LIBINT_SUPPORTS_ONEBODY
# * SUPPORT_T1G12

# "_candidate" variables are not needed for config.h but are used to figure out
#   the AM limits at the CMake level so that libint2-config.cmake components may
#   be defined and consuming codes can state their requirements. For example,
#   `find_package(Libint2 REQUIRED COMPONENTS eri_hhhh_d1)` requires the detected
#   library to include gradient integrals of at least AM=5.
#   See INSTALL.md for details.

# these bounds are for components. Above hard_max produces cmake warning. 0 produces
#   error in keeping with configure.ac logic.
set(LIBINT_HARD_MAX_AM 12)
set(LIBINT_HARD_MIN_AM 0)  # formerly 2
# amstr = "SPDFGHIKLMNOQRTUVWXYZ"
set(_am0 "s")
set(_am1 "p")
set(_am2 "d")
set(_am3 "f")
set(_am4 "g")
set(_am5 "h")
set(_am6 "i")
set(_am7 "k")
set(_am8 "l")
set(_am9 "m")
set(_am10 "n")
set(_am11 "o")
set(_am12 "q")
set(_AM0 "S")
set(_AM1 "P")
set(_AM2 "D")
set(_AM3 "F")
set(_AM4 "G")
set(_AM5 "H")
set(_AM6 "I")
set(_AM7 "K")
set(_AM8 "L")
set(_AM9 "M")
set(_AM10 "N")
set(_AM11 "O")
set(_AM12 "Q")

macro(numerical_max_of_list ansvar liste)
    set(_max "-100")
    foreach(_i ${liste})
        if (${_i} GREATER _max)
            set(_max "${_i}")
        endif()
    endforeach()
    set(${ansvar} "${_max}")
endmacro()


message(STATUS "Processing integrals classes ...")

# <<<  overall derivatives level  >>>

set(_glob_classes_derivs ${ENABLE_ONEBODY};${ENABLE_ERI};${ENABLE_ERI3};${ENABLE_ERI2};${ENABLE_G12};${ENABLE_G12DKH})
numerical_max_of_list(_max_deriv "${_glob_classes_derivs}")
message(STATUS "Preparing highest derivative level ${_max_deriv}")

# <<<  overall max_am defaults  >>>

list(LENGTH WITH_MAX_AM _ntokens_maxam)
if (_ntokens_maxam GREATER 1)
    math(EXPR _ntokens_xptd_max_deriv "${_max_deriv} + 1")
    if (NOT _ntokens_xptd_max_deriv EQUAL _ntokens_maxam)
        message(FATAL_ERROR "Invalid value for WITH_MAX_AM (${WITH_MAX_AM}). Highest ENABLE_ derivative (${_max_deriv}) requires list length ${_ntokens_xptd_max_deriv}, not ${_ntokens_maxam}.")
    endif()

    numerical_max_of_list(_max_am "${WITH_MAX_AM}")
    list(JOIN WITH_MAX_AM "," _sam)
    set(LIBINT_MAX_AM_LIST ${_sam})
    # when LIST populated, only overall LIBINT (not specific ints classes) sets both MAX_AM & MAX_AM_LIST
    set(LIBINT_MAX_AM ${_max_am})
    set(_max_LIBINT_MAX_AM ${LIBINT_MAX_AM})
else()
    set(LIBINT_MAX_AM_LIST "")
    set(LIBINT_MAX_AM ${WITH_MAX_AM})
endif()

foreach(_d RANGE 0 ${_max_deriv})
    if (${_d} LESS _ntokens_maxam)
        list(GET WITH_MAX_AM ${_d} _eri3_candidate0_d${_d})
        set(_dflt_candidate0_d${_d} "${LIBINT_MAX_AM}")
    else()
        set(_eri3_candidate0_d${_d} "${LIBINT_MAX_AM}")
        set(_dflt_candidate0_d${_d} "-1")
    endif()
    # _candidate0_dD=int_am defined up to highest ENABLE_cls deriv D from best info from WITH_MAX_AM
    message(VERBOSE "setting _eri3_candidate0_d${_d}=${_eri3_candidate0_d${_d}}")
    message(VERBOSE "setting _dflt_candidate0_d${_d}=${_dflt_candidate0_d${_d}}")
endforeach()

if (LIBINT_MAX_AM GREATER_EQUAL ${LIBINT_HARD_MAX_AM})
    message(WARNING "LIBINT_MAX_AM=${LIBINT_MAX_AM} is greater than ${LIBINT_HARD_MAX_AM}. Are you sure you know what you are doing?")
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
    numerical_max_of_list(LIBINT_OPT_AM "${_processed_OPT_AM_LIST}")
else()
    if(WITH_OPT_AM EQUAL -1)
        # first branch is a nice default pattern but not exactly what configure.ac prescribes, so bypassing it
        # if (_ntokens_maxam GREATER 1)
        if (FALSE)
            # no opt and list max: default list opt from max
            set(_processed_OPT_AM_LIST )
            math(EXPR _range_limit "${_ntokens_maxam} - 1")
            foreach(_d RANGE ${_range_limit})
                list(GET WITH_MAX_AM ${_d} _max_am)
                math(EXPR _opt_am "${_max_am}/2 + 1")
                list(APPEND _processed_OPT_AM_LIST ${_opt_am})
            endforeach()

            list(JOIN _processed_OPT_AM_LIST "," LIBINT_OPT_AM_LIST)
            numerical_max_of_list(LIBINT_OPT_AM "${_processed_OPT_AM_LIST}")
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
                set(_candidate0_${class}_d${_d} ${_dflt_candidate0_d${_d}})
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

                if (_candidate_${class}_d${_d} LESS_EQUAL 0)
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

                if (${class}_MAX_AM GREATER_EQUAL ${LIBINT_HARD_MAX_AM})
                    message(WARNING "Value for ${class}_MAX_AM too high (${${class}_MAX_AM} >= ${LIBINT_HARD_MAX_AM}). Are you sure you know what you are doing?")
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


macro(process_integrals_class_alt class)

    list(LENGTH ENABLE_${class} _ntokens)
    if (NOT _ntokens EQUAL 1)
        message(FATAL_ERROR "Invalid value for ENABLE_${class} (${ENABLE_${class}}). Use scalar of maximum derivative level, not list.")
    endif()

    if (ENABLE_${class} GREATER_EQUAL 0)
        set(INCLUDE_${class} ${ENABLE_${class}})

        foreach(_d RANGE 0 ${_max_deriv})
            if (${_d} LESS_EQUAL ${INCLUDE_${class}})
                # no per-d defaults. use energy
                set(_candidate0_${class}_d${_d} ${_dflt_candidate0_d0})
                message(VERBOSE "setting _candidate0_${class}_d${_d}=${_candidate0_${class}_d${_d}}")
            endif()
        endforeach()

        set(LIBINT_SUPPORTS_${class} yes)
        set(LIBINT_${class}_DERIV ${INCLUDE_${class}})
        message(STATUS "Enabling integrals class ${class} to derivative ${INCLUDE_${class}}")
    else()
        set(INCLUDE_${class} "-1")
        set(${class}_MAX_AM "")
        message(STATUS "Disabling integrals class ${class}")
    endif()

    if (ENABLE_${class} GREATER_EQUAL 0)
        list(LENGTH WITH_${class}_MAX_AM _ntokens)
        if (_ntokens GREATER 1)
            message(FATAL_ERROR "Invalid value for WITH_${class}_MAX_AM (${WITH_${class}_MAX_AM}). ENABLE_${class} derivative supports only scalar, not list length ${_ntokens}.")

        else()
            if (WITH_${class}_MAX_AM EQUAL -1)
                foreach(_d RANGE ${INCLUDE_${class}})
                    set(_candidate_${class}_d${_d} ${_candidate0_${class}_d0})
                    message(VERBOSE "setting _candidate_${class}_d${_d}=${_candidate_${class}_d0}")
                endforeach()

                set(_${class}_MAX_AM_pre "")
                set(${class}_MAX_AM ${_candidate0_${class}_d0})
                # note: unlike usual classes, C++ code seems to want class_MAX_AM set explicitly to config.h
            else()
                set(_${class}_MAX_AM_pre ${WITH_${class}_MAX_AM})
                set(${class}_MAX_AM ${WITH_${class}_MAX_AM})

                foreach(_d RANGE ${INCLUDE_${class}})
                    set(_candidate_${class}_d${_d} ${${class}_MAX_AM})
                    message(VERBOSE "setting _candidate_${class}_d${_d}=${_candidate_${class}_d${_d}}")
                endforeach()

                if (${class}_MAX_AM GREATER_EQUAL ${LIBINT_HARD_MAX_AM})
                    message(WARNING "Value for ${class}_MAX_AM too high (${${class}_MAX_AM} >= ${LIBINT_HARD_MAX_AM}). Are you sure you know what you are doing?")
                elseif (${class}_MAX_AM LESS_EQUAL 0)
                    message(FATAL_ERROR "Invalid value for ${class}_MAX_AM (${${class}_MAX_AM} <= 0).")
                endif()
            endif()
        endif()
        message(STATUS "Enabling integrals class ${class} to max AM ${_${class}_MAX_AM_pre} (else ${LIBINT_MAX_AM})")

        list(LENGTH WITH_${class}_OPT_AM _ntokens)
        if (_ntokens GREATER 1)
            message(FATAL_ERROR "Invalid value for WITH_${class}_OPT_AM (${WITH_${class}_OPT_AM}). ENABLE_${class} derivative supports only scalar, not list length ${_ntokens}.")

        else()
            if (WITH_${class}_OPT_AM EQUAL -1)
                set(_${class}_OPT_AM_pre "")
                set(${class}_OPT_AM ${LIBINT_OPT_AM})
                # note: unlike usual classes, C++ code seems to want class_MAX_AM set explicitly
            else()
                set(_${class}_OPT_AM_pre ${WITH_${class}_OPT_AM})
                set(${class}_OPT_AM ${WITH_${class}_OPT_AM})
            endif()
        endif()
        message(STATUS "Enabling integrals class ${class} to opt AM ${_${class}_OPT_AM_pre} (else ${LIBINT_OPT_AM})")
    endif()
endmacro()


process_integrals_class(ONEBODY)
process_integrals_class(ERI)
process_integrals_class(ERI3)
process_integrals_class(ERI2)
# unlike above, these classes (1) don't do AM_LIST and (2) require value in config.h if enabled
process_integrals_class_alt(G12)
process_integrals_class_alt(G12DKH)

if (ENABLE_G12 GREATER_EQUAL 0)
    set(SUPPORT_T1G12 ${ENABLE_T1G12_SUPPORT})
else()
    set(SUPPORT_T1G12 OFF)
endif()

add_feature_info(
  "general integral"
  "ON"
  "config.h: LIBINT_MAX_AM=${LIBINT_MAX_AM} LIBINT_MAX_AM_LIST=${LIBINT_MAX_AM_LIST} LIBINT_OPT_AM=${LIBINT_OPT_AM} LIBINT_OPT_AM_LIST=${LIBINT_OPT_AM_LIST}"
  )

# form list of active class + deriv + max_am strings to use in libint2-config.cmake
# * this generates components in same order as export/cmake/configuration-gen.py
set(Libint2_ERI_COMPONENTS "")
set(_amlist "")
set(_eri3_impure_sh "")
set(_eri2_impure_sh "")
add_feature_info(
  "integral class MULTIPOLE derivative 0"
  "MULTIPOLE_MAX_ORDER"
  "max_am ${MULTIPOLE_MAX_ORDER}"
  )
foreach(_l RANGE ${LIBINT_HARD_MIN_AM} ${MULTIPOLE_MAX_ORDER})
    # RANGE with count-down works but docs say behavior is undefined, so we count-up and reverse
    # RANGE starting at 2 avoids enumerating s, p
    set(_lbl "multipole_${_am${_l}}${_am${_l}}_d0")
    list(APPEND _amlist "${_lbl}")
endforeach()
list(REVERSE _amlist)
list(APPEND Libint2_ERI_COMPONENTS "${_amlist}")
message(VERBOSE "setting components ${_amlist}")

foreach(_cls ONEBODY;ERI;ERI3;ERI2;G12;G12DKH)
    if((_cls STREQUAL G12) OR (_cls STREQUAL G12DKH))
        add_feature_info(
          "integral class ${_cls}"
          "INCLUDE_${_cls} GREATER -1"
          "config.h: INCLUDE_${_cls}=${INCLUDE_${_cls}} ${_cls}_MAX_AM=${${_cls}_MAX_AM} ${_cls}_OPT_AM=${${_cls}_OPT_AM}"
          )
    else()
        add_feature_info(
          "integral class ${_cls}"
          "INCLUDE_${_cls} GREATER -1"
          "config.h: INCLUDE_${_cls}=${INCLUDE_${_cls}} ${_cls}_MAX_AM=${${_cls}_MAX_AM} ${_cls}_MAX_AM_LIST=${${_cls}_MAX_AM_LIST} ${_cls}_OPT_AM=${${_cls}_OPT_AM} ${_cls}_OPT_AM_LIST=${${_cls}_OPT_AM_LIST}"
          )
    endif()
    if(_cls STREQUAL "G12DKH")
        # add G12DKH below when it's granted components
        continue()
    endif()
    if (INCLUDE_${_cls} GREATER -1)
        foreach (_d RANGE 0 ${INCLUDE_${_cls}})
            add_feature_info(
              "integral class ${_cls} derivative ${_d}"
              "INCLUDE_${_cls} GREATER -1"
              "max_am ${_candidate_${_cls}_d${_d}}"
              )
            set(_amlist "")
            set(_pureamlist "")
            foreach(_l RANGE ${LIBINT_HARD_MIN_AM} ${_candidate_${_cls}_d${_d}})  # LIBINT_cls_MAX_AM[_LIST]
                if (_cls STREQUAL "ERI")
                    list(APPEND _amlist             "eri_${_am${_l}}${_am${_l}}${_am${_l}}${_am${_l}}_d${_d}")
                elseif (_cls STREQUAL "ERI2")
                    list(APPEND _amlist             "eri_${_AM${_l}}${_AM${_l}}_d${_d}")
                    list(APPEND _pureamlist         "eri_${_am${_l}}${_am${_l}}_d${_d}")
                elseif (_cls STREQUAL "ONEBODY")
                    list(APPEND _amlist         "onebody_${_am${_l}}${_am${_l}}_d${_d}")
                elseif (_cls STREQUAL "G12")
                    list(APPEND _amlist             "g12_${_am${_l}}${_am${_l}}${_am${_l}}${_am${_l}}_d${_d}")
                endif()
            endforeach()
            if (_cls STREQUAL "ERI3")
                foreach(_lfit RANGE ${LIBINT_HARD_MIN_AM} ${_candidate_${_cls}_d${_d}})  # LIBINT_ERI3_MAX_AM[_LIST], fitting
                    foreach(_lpr RANGE ${LIBINT_HARD_MIN_AM} ${_eri3_candidate0_d${_d}})  # LIBINT_MAX_AM[_LIST], paired
                        if (_lfit GREATER_EQUAL _lpr)
                            list(APPEND _amlist     "eri_${_am${_lpr}}${_am${_lpr}}${_AM${_lfit}}_d${_d}")
                            list(APPEND _pureamlist "eri_${_am${_lpr}}${_am${_lpr}}${_am${_lfit}}_d${_d}")
                        endif()
                    endforeach()
                endforeach()
            endif()
            list(REVERSE _amlist)
            list(APPEND Libint2_ERI_COMPONENTS "${_amlist}")
            message(VERBOSE "setting components ${_amlist}")
            list(REVERSE _pureamlist)
            if (_cls STREQUAL "ERI2")
                list(APPEND _eri2_impure_sh "${_pureamlist}")
            elseif (_cls STREQUAL "ERI3")
                list(APPEND _eri3_impure_sh "${_pureamlist}")
            endif()
        endforeach()
        if ((_cls STREQUAL "ERI3") AND NOT ERI3_PURE_SH)
            list(APPEND Libint2_ERI_COMPONENTS "${_eri3_impure_sh}")
            message(VERBOSE "setting components ${_eri3_impure_sh}")
        elseif ((_cls STREQUAL "ERI2") AND NOT ERI2_PURE_SH)
            list(APPEND Libint2_ERI_COMPONENTS "${_eri2_impure_sh}")
            message(VERBOSE "setting components ${_eri2_impure_sh}")
        endif()
    endif()
endforeach()
message(STATUS "Library will satisfy ERI AM components: ${Libint2_ERI_COMPONENTS}")
