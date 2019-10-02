# records whether user set option so the defaulting logic works

if (DEFINED LIBINT_MAX_AM)
    set(_user_LIBINT_MAX_AM 1)
endif()
if (DEFINED LIBINT_OPT_AM)
    set(_user_LIBINT_OPT_AM 1)
endif()

if (DEFINED ONEBODY_MAX_AM)
    set(_user_ONEBODY_MAX_AM 1)
endif()
if (DEFINED ONEBODY_OPT_AM)
    set(_user_ONEBODY_OPT_AM 1)
endif()

if (DEFINED ERI_MAX_AM)
    set(_user_ERI_MAX_AM 1)
endif()
if (DEFINED ERI_OPT_AM)
    set(_user_ERI_OPT_AM 1)
endif()

if (DEFINED ERI3_MAX_AM)
    set(_user_ERI3_MAX_AM 1)
endif()
if (DEFINED ERI3_OPT_AM)
    set(_user_ERI3_OPT_AM 1)
endif()

if (DEFINED ERI2_MAX_AM)
    set(_user_ERI2_MAX_AM 1)
endif()
if (DEFINED ERI2_OPT_AM)
    set(_user_ERI2_OPT_AM 1)
endif()


if (DEFINED G12_MAX_AM)
    set(_user_G12_MAX_AM 1)
endif()
if (DEFINED G12_OPT_AM)
    set(_user_G12_OPT_AM 1)
endif()


if (DEFINED G12DKH_MAX_AM)
    set(_user_G12DKH_MAX_AM 1)
endif()
if (DEFINED G12DKH_OPT_AM)
    set(_user_G12DKH_OPT_AM 1)
endif()
