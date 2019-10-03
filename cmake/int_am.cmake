# handle the defaulting and setting of the following variables
# * ENABLE_[ONEBODY|ERI|ERI2|ERI3|G12|G12DKH]
# * [LIBINT|ONEBODY|ERI|ERI2|ERI3|G12|G12DKH]_[MAX|OPT]_AM[|_LIST]
# * LIBINT_ONEBODY_DERIV
# * LIBINT_SUPPORTS_ONEBODY


message(STATUS "Processing integrals classes ...")

# <<<  overall defaults (LIBINT_MAX/OPT_AM)  >>>

list(LENGTH LIBINT_MAX_AM _lam)
if (_lam GREATER 1)
    list(JOIN LIBINT_MAX_AM "," _sam)
    execute_process (COMMAND bash -c "echo ${_sam} | tr , '\n' | sort -n | tail -n1"
                     OUTPUT_VARIABLE _max_LIBINT_MAX_AM)
    set(LIBINT_MAX_AM_LIST ${_sam})
    set(LIBINT_MAX_AM "")
else()
    set(_max_LIBINT_MAX_AM ${LIBINT_MAX_AM})

    if (LIBINT_MAX_AM GREATER_EQUAL 8)
        message(FATAL "Value for LIBINT_MAX_AM too high (${LIBINT_MAX_AM}). Are you sure you know what you are doing?")
    elseif (LIBINT_MAX_AM LESS_EQUAL 0)
        message(FATAL "Invalid value for LIBINT_MAX_AM (${LIBINT_MAX_AM}).")
    endif()
endif()

list(LENGTH LIBINT_OPT_AM _lam)
if (_lam GREATER 1)
    list(JOIN LIBINT_OPT_AM "," _sam)
    execute_process (COMMAND bash -c "echo ${_sam} | tr , '\n' | sort -n | tail -n1"
                     OUTPUT_VARIABLE _max_LIBINT_OPT_AM)
    set(LIBINT_OPT_AM_LIST ${_sam})
    set(LIBINT_OPT_AM "")
else()
    if (NOT _user_LIBINT_OPT_AM)
        math(EXPR LIBINT_OPT_AM "${_max_LIBINT_MAX_AM}/2 + 1")
        #math(EXPR LIBINT_OPT_AM "${_max_LIBINT_MAX_AM}/2 + ${_max_LIBINT_MAX_AM}%2")
    endif()

    if (LIBINT_OPT_AM GREATER LIBINT_MAX_AM)
        message(FATAL "Invalid value for LIBINT_OPT_AM (${LIBINT_OPT_AM} !<= ${LIBINT_MAX_AM}).")
    endif()
endif()


# <<<  Macro  >>>

macro(process_integrals_class class)
    if (ENABLE_${class} GREATER_EQUAL 0)
        set(INCLUDE_${class} ${ENABLE_${class}})
        if (${class} STREQUAL ONEBODY)
            set(LIBINT_SUPPORTS_${class} yes)
            set(LIBINT_${class}_DERIV ${INCLUDE_${class}})
        endif()
        message(STATUS "Enabling integrals class ${class} to derivative ${INCLUDE_${class}}")
    else()
        set(INCLUDE_${class} "-1")
        set(${class}_MAX_AM "")
        set(${class}_MAX_AM_LIST "")
        message(STATUS "Disabling integrals class ${class}")
    endif()
    
    if (ENABLE_${class} GREATER_EQUAL 0)
        list(LENGTH ${class}_MAX_AM _lam)
        if (_lam GREATER 1)
            list(JOIN ${class}_MAX_AM "," _sam)
            execute_process (COMMAND bash -c "echo ${_sam} | tr , '\n' | sort -n | tail -n1"
                             OUTPUT_VARIABLE _max_${class}_MAX_AM)
            set(${class}_MAX_AM_LIST ${_sam})
            set(${class}_MAX_AM "")
        else()
            if (NOT _user_${class}_MAX_AM)
                set(${class}_MAX_AM "")
                set(_max_${class}_MAX_AM ${LIBINT_MAX_AM})
            else()
                # _max_* variable in case want to default opt_am from some day
                set(_max_${class}_MAX_AM ${${class}_MAX_AM})
                if (${class}_MAX_AM GREATER_EQUAL 8)
                    message(FATAL "Value for ${class}_MAX_AM too high (${${class}_MAX_AM}). Are you sure you know what you are doing?")
                elseif (${class}_MAX_AM LESS_EQUAL 0)
                    message(FATAL "Invalid value for ${class}_MAX_AM (${${class}_MAX_AM}).")
                endif()
            endif()
        endif()
        message(STATUS "Enabling integrals class ${class} to max AM ${${class}_MAX_AM}${${class}_MAX_AM_LIST} (else ${LIBINT_MAX_AM})")
    
        list(LENGTH ${class}_OPT_AM _lam)
        if (_lam GREATER 1)
            list(JOIN ${class}_OPT_AM "," _sam)
            execute_process (COMMAND bash -c "echo ${_sam} | tr , '\n' | sort -n | tail -n1"
                             OUTPUT_VARIABLE _max_${class}_OPT_AM)
            set(${class}_OPT_AM_LIST ${_sam})
            set(${class}_OPT_AM "")
        else()
            if (_user_${class}_OPT_AM)
                if (${class}_OPT_AM GREATER ${class}_MAX_AM)
                    message(FATAL "Invalid value for ${class}_OPT_AM (${${class}_OPT_AM} !<= ${${class}_MAX_AM}).")
                endif()
            else()
                set(${class}_OPT_AM "")
            endif()
    
        endif()
        message(STATUS "Enabling integrals class ${class} to opt AM ${${class}_OPT_AM}${${class}_OPT_AM_LIST} (else ${LIBINT_OPT_AM})")
    endif()
endmacro()


process_integrals_class(ONEBODY)
process_integrals_class(ERI)
process_integrals_class(ERI3)
process_integrals_class(ERI2)

# discrepancy, as configure doesn't do AM_LIST for these
process_integrals_class(G12)
process_integrals_class(G12DKH)
