
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


# <<<  ERI  >>>

if (ENABLE_ERI GREATER_EQUAL 0)
    set(INCLUDE_ERI ${ENABLE_ERI})
else()
    set(INCLUDE_ERI "-1")
endif()
message("ERI ${INCLUDE_ERI}")

list(LENGTH ERI_MAX_AM _lam)
#message("ERI_MAX_AM ${ERI_MAX_AM} ${_lam}")
if (_lam GREATER 1)
    list(JOIN ERI_MAX_AM "," _sam)
    #message("_sam ${_sam}")
    execute_process (COMMAND bash -c "echo ${_sam} | tr , '\n' | sort -n | tail -n1"
                     OUTPUT_VARIABLE _max_ERI_MAX_AM)
    set(ERI_MAX_AM_LIST ${_sam})
    set(ERI_MAX_AM "")
    #message("LIST ${ERI_MAX_AM} ${ERI_MAX_AM_LIST}")
else()
    if (ERI_MAX_AM GREATER_EQUAL 8)
        message(FATAL "Value for ERI_MAX_AM too high (${ERI_MAX_AM}). Are you sure you know what you are doing?")
    elseif (ERI_MAX_AM LESS_EQUAL 0)
        message(FATAL "Invalid value for ERI_MAX_AM (${ERI_MAX_AM}).")
    endif()
    #message("SCAL ${ERI_MAX_AM} ${ERI_MAX_AM_LIST}")
    set(_max_ERI_MAX_AM ${ERI_MAX_AM})
endif()

list(LENGTH ERI_OPT_AM _lam)
#message("ERI_OPT_AM ${ERI_OPT_AM} ${_lam}")
if (_lam GREATER 1)
    list(JOIN ERI_OPT_AM "," _sam)
    #message("_sam ${_sam}")
    execute_process (COMMAND bash -c "echo ${_sam} | tr , '\n' | sort -n | tail -n1"
                     OUTPUT_VARIABLE _max_ERI_OPT_AM)
    set(ERI_OPT_AM_LIST ${_sam})
    set(ERI_OPT_AM "")
    #message("LIST ${ERI_OPT_AM} ${ERI_OPT_AM_LIST}")
else()
    if (ERI_OPT_AM EQUAL -1)
        #math(EXPR ERI_OPT_AM "${_max_ERI_MAX_AM}/2 + ${_max_ERI_MAX_AM}%2")
        math(EXPR ERI_OPT_AM "${_max_ERI_MAX_AM}/2 + 1")
    endif()
    if (ERI_OPT_AM GREATER ERI_MAX_AM)
        message(FATAL "Invalid value for ERI_OPT_AM (${ERI_OPT_AM} !<= ${ERI_MAX_AM}).")
    endif()
    #message("SCAL ${ERI_OPT_AM} ${ERI_OPT_AM_LIST}")
endif()

# temp
set(LIBINT_MAX_AM ${_max_ERI_MAX_AM})
set(LIBINT_OPT_AM ${ERI_OPT_AM})

# <<<  ERI3  >>>

if (ENABLE_ERI3 GREATER_EQUAL 0)
    set(INCLUDE_ERI3 ${ENABLE_ERI3})
else()
    set(INCLUDE_ERI3 "-1")
endif()
message("ERI3 ${INCLUDE_ERI3}")


# <<<  ERI2  >>>

if (ENABLE_ERI2 GREATER_EQUAL 0)
    set(INCLUDE_ERI2 ${ENABLE_ERI2})
else()
    set(INCLUDE_ERI2 "-1")
endif()
message("ERI2 ${INCLUDE_ERI2}")
