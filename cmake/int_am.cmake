
#LIBINT_OPT_AM_LIST
#LIBINT_MAX_AM_LIST
#
#ERI_OPT_AM
#ERI_OPT_AM_LIST
#ERI_MAX_AM
#ERI_MAX_AM_LIST
#
#ERI2_OPT_AM
#ERI2_OPT_AM_LIST
#ERI2_MAX_AM
#ERI2_MAX_AM_LIST
#
#ERI3_OPT_AM
#ERI3_OPT_AM_LIST
#ERI3_MAX_AM
#ERI3_MAX_AM_LIST
#
#LIBINT_ONEBODY_DERIV
#LIBINT_SUPPORTS_ONEBODY
#ONEBODY_OPT_AM
#ONEBODY_OPT_AM_LIST
#ONEBODY_MAX_AM
#ONEBODY_MAX_AM_LIST
#
#G12_OPT_AM
#G12_MAX_AM
#
#G12DKH_OPT_AM
#G12DKH_MAX_AM


# <<<  LIBINT  >>>

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
#message(STATUS "LIBINT integrals enabled to max AM ${LIBINT_MAX_AM} ${LIBINT_MAX_AM_LIST}")

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
    endif()

    if (LIBINT_OPT_AM GREATER LIBINT_MAX_AM)
        message(FATAL "Invalid value for LIBINT_OPT_AM (${LIBINT_OPT_AM} !<= ${LIBINT_MAX_AM}).")
    endif()
endif()
#message(STATUS "LIBINT integrals enabled to opt AM ${LIBINT_OPT_AM} ${LIBINT_OPT_AM_LIST}")


# <<<  ONEBODY  >>>

if (ENABLE_ONEBODY GREATER_EQUAL 0)
    set(INCLUDE_ONEBODY ${ENABLE_ONEBODY})
    set(LIBINT_SUPPORTS_ONEBODY yes)
    set(LIBINT_ONEBODY_DERIV ${INCLUDE_ONEBODY})
    message(STATUS "Enabling integrals class ONEBODY to derivative ${INCLUDE_ONEBODY}")
else()
    set(INCLUDE_ONEBODY "-1")
    set(ONEBODY_MAX_AM "")
    set(ONEBODY_MAX_AM_LIST "")
    message(STATUS "Disabling integrals class ONEBODY")
endif()

if (ENABLE_ONEBODY GREATER_EQUAL 0)
    list(LENGTH ONEBODY_MAX_AM _lam)
    if (_lam GREATER 1)
        list(JOIN ONEBODY_MAX_AM "," _sam)
        execute_process (COMMAND bash -c "echo ${_sam} | tr , '\n' | sort -n | tail -n1"
                         OUTPUT_VARIABLE _max_ONEBODY_MAX_AM)
        set(ONEBODY_MAX_AM_LIST ${_sam})
        set(ONEBODY_MAX_AM "")
    else()
        if (NOT _user_ONEBODY_MAX_AM)
            set(ONEBODY_MAX_AM "")
            set(_max_ONEBODY_MAX_AM ${LIBINT_MAX_AM})
        else()
            set(_max_ONEBODY_MAX_AM ${ONEBODY_MAX_AM})
            if (ONEBODY_MAX_AM GREATER_EQUAL 8)
                message(FATAL "Value for ONEBODY_MAX_AM too high (${ONEBODY_MAX_AM}). Are you sure you know what you are doing?")
            elseif (ONEBODY_MAX_AM LESS_EQUAL 0)
                message(FATAL "Invalid value for ONEBODY_MAX_AM (${ONEBODY_MAX_AM}).")
            endif()
        endif()
    endif()
    message(STATUS "Enabling integrals class ONEBODY to max AM ${ONEBODY_MAX_AM}${ONEBODY_MAX_AM_LIST}")

    list(LENGTH ONEBODY_OPT_AM _lam)
    if (_lam GREATER 1)
        list(JOIN ONEBODY_OPT_AM "," _sam)
        execute_process (COMMAND bash -c "echo ${_sam} | tr , '\n' | sort -n | tail -n1"
                         OUTPUT_VARIABLE _max_ONEBODY_OPT_AM)
        set(ONEBODY_OPT_AM_LIST ${_sam})
        set(ONEBODY_OPT_AM "")
    else()
        if (_user_ONEBODY_OPT_AM)
            if (ONEBODY_OPT_AM GREATER ONEBODY_MAX_AM)
                message(FATAL "Invalid value for ONEBODY_OPT_AM (${ONEBODY_OPT_AM} !<= ${ONEBODY_MAX_AM}).")
            endif()
        else()
            set(ONEBODY_OPT_AM "")
        endif()
    #        #math(EXPR ONEBODY_OPT_AM "${_max_ONEBODY_MAX_AM}/2 + ${_max_ONEBODY_MAX_AM}%2")
    #        math(EXPR ONEBODY_OPT_AM "${_max_ONEBODY_MAX_AM}/2 + 1")

    endif()
    message(STATUS "Enabling integrals class ONEBODY to opt AM ${ONEBODY_OPT_AM}${ONEBODY_OPT_AM_LIST}")
endif()


# <<<  ERI  >>>

if (ENABLE_ERI GREATER_EQUAL 0)
    set(INCLUDE_ERI ${ENABLE_ERI})
    message(STATUS "Enabling integrals class ERI to derivative ${INCLUDE_ERI}")
else()
    set(INCLUDE_ERI "-1")
    set(ERI_MAX_AM "")
    set(ERI_MAX_AM_LIST "")
    message(STATUS "Disabling integrals class ERI")
endif()

if (ENABLE_ERI GREATER_EQUAL 0)
    list(LENGTH ERI_MAX_AM _lam)
    if (_lam GREATER 1)
        list(JOIN ERI_MAX_AM "," _sam)
        execute_process (COMMAND bash -c "echo ${_sam} | tr , '\n' | sort -n | tail -n1"
                         OUTPUT_VARIABLE _max_ERI_MAX_AM)
        set(ERI_MAX_AM_LIST ${_sam})
        set(ERI_MAX_AM "")
    else()
        if (NOT _user_ERI_MAX_AM)
            set(ERI_MAX_AM ${LIBINT_MAX_AM})
        endif()
        set(_max_ERI_MAX_AM ${ERI_MAX_AM})

        if (ERI_MAX_AM GREATER_EQUAL 8)
            message(FATAL "Value for ERI_MAX_AM too high (${ERI_MAX_AM}). Are you sure you know what you are doing?")
        elseif (ERI_MAX_AM LESS_EQUAL 0)
            message(FATAL "Invalid value for ERI_MAX_AM (${ERI_MAX_AM}).")
        endif()
    endif()
    message(STATUS "Enabling integrals class ERI to max AM ${ERI_MAX_AM}${ERI_MAX_AM_LIST}")

    list(LENGTH ERI_OPT_AM _lam)
    if (_lam GREATER 1)
        list(JOIN ERI_OPT_AM "," _sam)
        execute_process (COMMAND bash -c "echo ${_sam} | tr , '\n' | sort -n | tail -n1"
                         OUTPUT_VARIABLE _max_ERI_OPT_AM)
        set(ERI_OPT_AM_LIST ${_sam})
        set(ERI_OPT_AM "")
    else()
        if (NOT _user_ERI_OPT_AM)
            set(ERI_OPT_AM ${LIBINT_OPT_AM})
        endif()
    #        #math(EXPR ERI_OPT_AM "${_max_ERI_MAX_AM}/2 + ${_max_ERI_MAX_AM}%2")
    #        math(EXPR ERI_OPT_AM "${_max_ERI_MAX_AM}/2 + 1")

        if (ERI_OPT_AM GREATER ERI_MAX_AM)
            message(FATAL "Invalid value for ERI_OPT_AM (${ERI_OPT_AM} !<= ${ERI_MAX_AM}).")
        endif()
    endif()
    message(STATUS "Enabling integrals class ERI to opt AM ${ERI_OPT_AM}${ERI_OPT_AM_LIST}")
endif()


# <<<  ERI3  >>>

if (ENABLE_ERI3 GREATER_EQUAL 0)
    set(INCLUDE_ERI3 ${ENABLE_ERI3})
    message(STATUS "Enabling integrals class ERI3 to derivative ${INCLUDE_ERI3}")
else()
    set(INCLUDE_ERI3 "-1")
    set(ERI3_MAX_AM "")
    set(ERI3_MAX_AM_LIST "")
    message(STATUS "Disabling integrals class ERI3")
endif()

if (ENABLE_ERI3 GREATER_EQUAL 0)
    list(LENGTH ERI3_MAX_AM _lam)
    if (_lam GREATER 1)
        list(JOIN ERI3_MAX_AM "," _sam)
        execute_process (COMMAND bash -c "echo ${_sam} | tr , '\n' | sort -n | tail -n1"
                         OUTPUT_VARIABLE _max_ERI3_MAX_AM)
        set(ERI3_MAX_AM_LIST ${_sam})
        set(ERI3_MAX_AM "")
    else()
        if (NOT _user_ERI3_MAX_AM)
            set(ERI3_MAX_AM ${LIBINT_MAX_AM})
        endif()
        set(_max_ERI3_MAX_AM ${ERI3_MAX_AM})

        if (ERI3_MAX_AM GREATER_EQUAL 8)
            message(FATAL "Value for ERI3_MAX_AM too high (${ERI3_MAX_AM}). Are you sure you know what you are doing?")
        elseif (ERI3_MAX_AM LESS_EQUAL 0)
            message(FATAL "Invalid value for ERI3_MAX_AM (${ERI3_MAX_AM}).")
        endif()
    endif()
    message(STATUS "Enabling integrals class ERI3 to max AM ${ERI3_MAX_AM}${ERI3_MAX_AM_LIST}")

    list(LENGTH ERI3_OPT_AM _lam)
    if (_lam GREATER 1)
        list(JOIN ERI3_OPT_AM "," _sam)
        execute_process (COMMAND bash -c "echo ${_sam} | tr , '\n' | sort -n | tail -n1"
                         OUTPUT_VARIABLE _max_ERI3_OPT_AM)
        set(ERI3_OPT_AM_LIST ${_sam})
        set(ERI3_OPT_AM "")
    else()
        if (NOT _user_ERI3_OPT_AM)
            set(ERI3_OPT_AM ${LIBINT_OPT_AM})
        endif()
    #        #math(EXPR ERI3_OPT_AM "${_max_ERI3_MAX_AM}/2 + ${_max_ERI3_MAX_AM}%2")
    #        math(EXPR ERI3_OPT_AM "${_max_ERI3_MAX_AM}/2 + 1")

        if (ERI3_OPT_AM GREATER ERI3_MAX_AM)
            message(FATAL "Invalid value for ERI3_OPT_AM (${ERI3_OPT_AM} !<= ${ERI3_MAX_AM}).")
        endif()
    endif()
    message(STATUS "Enabling integrals class ERI3 to opt AM ${ERI3_OPT_AM}${ERI3_OPT_AM_LIST}")
endif()


# <<<  ERI2  >>>

if (ENABLE_ERI2 GREATER_EQUAL 0)
    set(INCLUDE_ERI2 ${ENABLE_ERI2})
    message(STATUS "Enabling integrals class ERI2 to derivative ${INCLUDE_ERI2}")
else()
    set(INCLUDE_ERI2 "-1")
    set(ERI2_MAX_AM "")
    set(ERI2_MAX_AM_LIST "")
    message(STATUS "Disabling integrals class ERI2")
endif()


# <<<  G12  >>>

if (ENABLE_G12 GREATER_EQUAL 0)
    set(INCLUDE_G12 ${ENABLE_G12})
    message(STATUS "Enabling integrals class G12 to derivative ${INCLUDE_G12}")
else()
    set(INCLUDE_G12 "-1")
    set(G12_MAX_AM "")
    set(G12_MAX_AM_LIST "")
    message(STATUS "Disabling integrals class G12")
endif()


# <<<  G12DKH  >>>

if (ENABLE_G12DKH GREATER_EQUAL 0)
    set(INCLUDE_G12DKH ${ENABLE_G12DKH})
    message(STATUS "Enabling integrals class G12DKH to derivative ${INCLUDE_G12DKH}")
else()
    set(INCLUDE_G12DKH "-1")
    set(G12DKH_MAX_AM "")
    set(G12DKH_MAX_AM_LIST "")
    message(STATUS "Disabling integrals class G12DKH")
endif()
