# Keep libtool and cmake in sync by collecting versions from configure.ac
# If CMake becomes primary buildsys, define as `project(..., VERSION)`

# <<<  Build Version  >>>

file(STRINGS "configure.ac" _libint_configure_ac
     REGEX "libint_mmm_version")
foreach(ver ${_libint_configure_ac})
    if (ver MATCHES "^define..libint_mmm_version...([0-9]+).([0-9]+).([0-9]+)..$")
        set(LIBINT_MAJOR_VERSION ${CMAKE_MATCH_1})
        set(LIBINT_MINOR_VERSION ${CMAKE_MATCH_2})
        set(LIBINT_MICRO_VERSION ${CMAKE_MATCH_3})
    endif()
endforeach()

set(Libint2_VERSION ${LIBINT_MAJOR_VERSION}.${LIBINT_MINOR_VERSION}.${LIBINT_MICRO_VERSION})


# <<<  Dev Version  >>>

file(STRINGS "configure.ac" _libint_configure_ac
     REGEX "libint_buildid")
foreach(ver ${_libint_configure_ac})
    if (ver MATCHES "^define..libint_buildid...([a-z0-9.]+)..$")
        set(LIBINT_TWEAK_VERSION ${CMAKE_MATCH_1})
    endif()
endforeach()

set(LIBINT_VERSION ${Libint2_VERSION}-${LIBINT_TWEAK_VERSION})
message(STATUS "Libint Full ${LIBINT_VERSION} Numeric ${Libint2_VERSION}")


# <<<  ABI Version  >>>

file(STRINGS "configure.ac" _libint_configure_ac
     REGEX "libint_so_version")
foreach(ver ${_libint_configure_ac})
    if (ver MATCHES "^define..libint_so_version...([0-9]+).([0-9]+).([0-9]+)..$")
        set(LIBINT_CURRENT_SOVERSION ${CMAKE_MATCH_1})
        set(LIBINT_REVISION_SOVERSION ${CMAKE_MATCH_2})
        set(LIBINT_AGE_SOVERSION ${CMAKE_MATCH_3})
    endif()
endforeach()

set(LIBINT_SOVERSION ${LIBINT_CURRENT_SOVERSION}:${LIBINT_REVISION_SOVERSION}:${LIBINT_AGE_SOVERSION})
math(EXPR LIBINT_MAJOR_SOVERSION "${LIBINT_CURRENT_SOVERSION} - ${LIBINT_AGE_SOVERSION}")
message(STATUS "Libint SO Full ${LIBINT_SOVERSION} Major ${LIBINT_MAJOR_SOVERSION}")
