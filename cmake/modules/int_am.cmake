# handle the defaulting and setting of the following variables
# * ENABLE_[ONEBODY|ERI2|ERI3|ERI|G12|G12DKH]
# * [LIBINT|ONEBODY|ERI2|ERI3|ERI|G12|G12DKH]_[MAX|OPT]_AM[|_LIST]
# * LIBINT_ONEBODY_DERIV
# * LIBINT_SUPPORTS_ONEBODY

# _candidate variables not needed for config.h but are used to figure
#   out the AM levels at the CMake level
#   so that libint2-config components may be defined and client codes can
#   require the detected library include gradient integrals up to AM=5 with
#   `find_package(Libint2 COMPONENTS g5)`


message(STATUS "Processing integrals classes ...")

# <<<  overall defaults (LIBINT_MAX/OPT_AM)  >>>

message(STATUS "WITH_MAX_AM=${WITH_MAX_AM}")
list(LENGTH WITH_MAX_AM _ntokens)
if (_ntokens GREATER 1)

    set(_given_max_am_list TRUE)
    set(_max_am 0)
    math(EXPR _max_deriv "${_ntokens} - 1")
    foreach(_d RANGE 0 ${_max_deriv})
        list(GET WITH_MAX_AM ${_d} _candidate0_d${_d})
        if (${_candidate0_d${_d}} GREATER _max_am)
            set(_max_am ${_candidate0_d${_d}})
        endif()
    endforeach()

    list(JOIN WITH_MAX_AM "," _sam)
    set(LIBINT_MAX_AM_LIST ${_sam})
    set(LIBINT_MAX_AM ${_max_am})  # only overall LIBINT, not specific integrals classes, sets both MAX_AM & MAX_AM_LIST
else()
    set(_given_max_am_list FALSE)
    set(_candidate0_d0 ${WITH_MAX_AM})

    set(LIBINT_MAX_AM_LIST "")
    set(LIBINT_MAX_AM ${WITH_MAX_AM})
    set(_max_LIBINT_MAX_AM ${LIBINT_MAX_AM})

    if (LIBINT_MAX_AM GREATER_EQUAL 8)
        message(FATAL "LIBINT_MAX_AM=${LIBINT_MAX_AM} is greater than 8. Are you sure you know what you are doing?")
    elseif (LIBINT_MAX_AM LESS_EQUAL 0)
        message(FATAL "Invalid value for LIBINT_MAX_AM (${LIBINT_MAX_AM}).")
    endif()
endif()

message(STATUS "LIBINT_MAX_AM_LIST=${LIBINT_MAX_AM_LIST} LIBINT_MAX_AM=${LIBINT_MAX_AM}")

list(LENGTH WITH_OPT_AM _ntokens)
if (_ntokens GREATER 1)
    set(_max_am 0)
    set(PROCESSED_OPT_AM_LIST )
    math(EXPR _max_deriv "${_ntokens} - 1")
    foreach(_d RANGE 0 ${_max_deriv})
        list(GET WITH_OPT_AM ${_d} _opt_am)
        if (_opt_am GREATER _max_opt_am)
            set(_max_opt_am ${_opt_am})
        endif()
        if (_have_max_am_list)
            list(GET WITH_MAX_AM ${_d} _max_am)
            if (_opt_am GREATER _max_am)
                list(APPEND PROCESSED_OPT_AM_LIST ${_max_am})
            else()
                list(APPEND PROCESSED_OPT_AM_LIST ${_opt_am})
            endif()
        else()
            if (_opt_am GREATER LIBINT_MAX_AM)
                list(APPEND PROCESSED_OPT_AM_LIST ${LIBINT_MAX_AM})
            else()
                list(APPEND PROCESSED_OPT_AM_LIST ${LIBINT_MAX_AM})
            endif()
        endif()

    endforeach()
    list(JOIN PROCESSED_OPT_AM_LIST "," _sam)
    set(LIBINT_OPT_AM_LIST ${_sam})
    set(LIBINT_OPT_AM ${_max_opt_am})  # only overall LIBINT, not specific integrals classes, sets both OPT_AM & OPT_AM_LIST
else()
    set(LIBINT_OPT_AM_LIST "")
    message(STATUS "WITH_OPT_AM=${WITH_OPT_AM}")
    if (WITH_OPT_AM EQUAL -1)
        math(EXPR LIBINT_OPT_AM "${_max_LIBINT_MAX_AM}/2 + 1")
        #math(EXPR LIBINT_OPT_AM "${_max_LIBINT_MAX_AM}/2 + ${_max_LIBINT_MAX_AM}%2")
    else()
        set(LIBINT_OPT_AM ${WITH_OPT_AM})
    endif()

    if (LIBINT_OPT_AM GREATER _max_LIBINT_MAX_AM)
        set(LIBINT_OPT_AM ${_max_LIBINT_MAX_AM})
    endif()
endif()

message(STATUS "LIBINT_OPT_AM_LIST=${LIBINT_OPT_AM_LIST} LIBINT_OPT_AM=${LIBINT_OPT_AM}")

# <<<  Macro  >>>

macro(process_integrals_class class)
    if (ENABLE_${class} GREATER_EQUAL 0)
        set(INCLUDE_${class} ${ENABLE_${class}})
        set(_candidate0_${class}_d${INCLUDE_${class}} ${_candidate0_d${INCLUDE_${class}}})
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
            math(EXPR _max_deriv "${_ntokens} - 1")
            foreach(_d RANGE 0 ${_max_deriv})
                list(GET WITH_${class}_MAX_AM _d _candidate_${class}_d${_d})
            endforeach()

            list(JOIN WITH_${class}_MAX_AM "," _sam)
            set(${class}_MAX_AM_LIST ${_sam})
            set(${class}_MAX_AM "")
        else()
            set(${class}_MAX_AM_LIST "")
            if (WITH_${class}_MAX_AM EQUAL -1)

                set(_candidate_${class}_E ${_candidate0_${class}_E})
                if (${INCLUDE_${class}} GREATER_EQUAL 1)
                    if (${_candidate0_${class}_G} EQUAL -1)
                        set(_candidate_${class}_G ${_candidate0_${class}_E})
                    else()
                        set(_candidate_${class}_G ${_candidate0_${class}_G})
                    endif()
                endif()
                if (${INCLUDE_${class}} GREATER_EQUAL 2)
                    if (${_candidate0_${class}_H} EQUAL -1)
                        set(_candidate_${class}_H ${_candidate0_${class}_E})
                    else()
                        set(_candidate_${class}_H ${_candidate0_${class}_H})
                    endif()
                endif()

                if (${INCLUDE_${class}} GREATER_EQUAL 0)
                    set(${class}_MAX_AM ${LIBINT_MAX_AM})
                else()
                    set(${class}_MAX_AM "")
                endif()
                set(_max_${class}_MAX_AM ${LIBINT_MAX_AM})
            else()
                # _max_* variable in case want to default opt_am from it some day
                set(${class}_MAX_AM ${WITH_${class}_MAX_AM})
                set(_max_${class}_MAX_AM ${${class}_MAX_AM})

                set(_candidate_${class}_E ${${class}_MAX_AM})
                if (${INCLUDE_${class}} GREATER_EQUAL 1)
                    set(_candidate_${class}_G ${${class}_MAX_AM})
                endif()
                if (${INCLUDE_${class}} GREATER_EQUAL 2)
                    set(_candidate_${class}_H ${${class}_MAX_AM})
                endif()

                if (${class}_MAX_AM GREATER_EQUAL 8)
                    message(FATAL "Value for ${class}_MAX_AM too high (${${class}_MAX_AM}). Are you sure you know what you are doing?")
                elseif (${class}_MAX_AM LESS_EQUAL 0)
                    message(FATAL "Invalid value for ${class}_MAX_AM (${${class}_MAX_AM}).")
                endif()
            endif()
        endif()
        message(STATUS "Enabling integrals class ${class} to max AM ${${class}_MAX_AM}${${class}_MAX_AM_LIST} (else ${LIBINT_MAX_AM}${LIBINT_MAX_AM_LIST})")

        list(LENGTH WITH_${class}_OPT_AM _lam)
        if (_lam GREATER 1)
            list(JOIN WITH_${class}_OPT_AM "," _sam)
            execute_process (COMMAND bash -c "echo ${_sam} | tr , '\n' | sort -n | tail -n1"
                             OUTPUT_VARIABLE _max_${class}_OPT_AM)
            set(${class}_OPT_AM_LIST ${_sam})
            set(${class}_OPT_AM "")
        else()
            set(${class}_OPT_AM_LIST "")
            if (WITH_${class}_OPT_AM EQUAL -1)
                set(${class}_OPT_AM "")
            else()
                set(${class}_OPT_AM ${WITH_${class}_OPT_AM})

                if (${class}_OPT_AM GREATER _max_${class}_MAX_AM)
                    message(FATAL "Invalid value for ${class}_OPT_AM (${${class}_OPT_AM} !<= ${_max_${class}_MAX_AM}).")
                endif()
            endif()

        endif()
        message(STATUS "Enabling integrals class ${class} to opt AM ${${class}_OPT_AM}${${class}_OPT_AM_LIST} (else ${LIBINT_OPT_AM}${LIBINT_OPT_AM_LIST})")
    endif()
endmacro()


process_integrals_class(ONEBODY)
process_integrals_class(ERI2)
process_integrals_class(ERI3)
process_integrals_class(ERI)

# discrepancy, as configure doesn't do AM_LIST for these
process_integrals_class(G12)
process_integrals_class(G12DKH)

# form list of active <class>_<deriv><max_am> strings to use in Libint2Config
set(Libint2_ERI_COMPONENTS "")
foreach(_cls ERI2;ERI3;ERI4)
    string(TOLOWER ${_cls} _lbl)
    set(_lbl "${_lbl}_")

    foreach (_deriv RANGE 0 ${INCLUDE_${_cls}})
        foreach(_l RANGE 0 ${_candidate_${_cls}_d${_deriv}})
            list(APPEND Libint2_ERI_COMPONENTS "${_lbl}_d${_deriv}_l${_l}")
        endforeach()
    endforeach()
endforeach()
message(STATUS "Library will satisfy ERI AM components: ${Libint2_ERI_COMPONENTS}")

