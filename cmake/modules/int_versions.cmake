# top-level CMakeLists.txt has defined:
# * PROJECT_VERSION_{MAJOR|MINOR|PATCH} through `project(... VERSION)`
# * Libint2Compiler_DESCRIPTION through `project(... DESCRIPTION)`
# * LIBINT_BUILDID
# * LIBINT_SO_VERSION
# * LIBINT_DOI
# * LibintRepository_{VERSION|DESCRIBE|COMMIT|DISTANCE} through `dynamic_version()`

# note that 3rd version integer is PATCH in CMake and MICRO in Libint
# see also int_computed.cmake.in for transmitting these values to the library's CMake
# note that with the dynamic/git scheme, it's important for repo to be up-to-date with tags


# <<<  Sortable Version  >>>

message(DEBUG "LibintRepository_VERSION  ${LibintRepository_VERSION}")
message(DEBUG "LibintRepository_COMMIT   ${LibintRepository_COMMIT}")
message(DEBUG "LibintRepository_DISTANCE ${LibintRepository_DISTANCE}")
message(DEBUG "LibintRepository_DESCRIBE ${LibintRepository_DESCRIBE}")

if (LibintRepository_DISTANCE STREQUAL "0")
    set(LIBINT_SORTABLE_VERSION "${LibintRepository_VERSION}")
else()
    set(LIBINT_SORTABLE_VERSION "${LibintRepository_VERSION}.post${LibintRepository_DISTANCE}")
endif()

string(SUBSTRING ${LibintRepository_COMMIT} 0 7 LIBINT_GIT_COMMIT)
message(DEBUG "LIBINT_GIT_COMMIT         ${LIBINT_GIT_COMMIT}")

# Below goes into BibTeX citation. Currently year of export. For year of tag, parse:
# `git show -s --no-notes --date=short --pretty='%cd' v2.7.2` responds: 2022-06-20
string(TIMESTAMP LIBINT_VERSION_YEAR "%Y")
message(DEBUG "LIBINT_VERSION_YEAR       ${LIBINT_VERSION_YEAR}")

set(LIBINT_DESCRIPTION "${Libint2Compiler_DESCRIPTION}")
message(DEBUG "LIBINT_DESCRIPTION        ${LIBINT_DESCRIPTION}")

# <<<  Build Version  >>>

set(LIBINT_MAJOR_VERSION ${PROJECT_VERSION_MAJOR})
set(LIBINT_MINOR_VERSION ${PROJECT_VERSION_MINOR})
set(LIBINT_MICRO_VERSION ${PROJECT_VERSION_PATCH})

set(LIBINT_VERSION ${LIBINT_MAJOR_VERSION}.${LIBINT_MINOR_VERSION}.${LIBINT_MICRO_VERSION})


# <<<  Dev Version  >>>

if (LIBINT_BUILDID)
  set(LIBINT_EXT_VERSION ${LIBINT_VERSION}-${LIBINT_BUILDID})
else()
  set(LIBINT_EXT_VERSION ${LIBINT_VERSION})
endif()

message(STATUS "Version: Full ${LIBINT_EXT_VERSION} Numeric ${LIBINT_VERSION} Sortable ${LIBINT_SORTABLE_VERSION}")

if (NOT(LibintRepository_VERSION STREQUAL LIBINT_VERSION))
    message(AUTHOR_WARNING
        "Version processing has gone wrong: ${LibintRepository_VERSION} != ${LIBINT_VERSION}")
endif()


# <<<  ABI Version  >>>

string(REPLACE ":" ";" LIBINT_SO_VERSION_LIST ${LIBINT_SO_VERSION})

list(GET LIBINT_SO_VERSION_LIST 0 LIBINT_CURRENT_SO_VERSION)
list(GET LIBINT_SO_VERSION_LIST 1 LIBINT_REVISION_SO_VERSION)
list(GET LIBINT_SO_VERSION_LIST 2 LIBINT_AGE_SO_VERSION)

math(EXPR LIBINT_MAJOR_SO_VERSION "${LIBINT_CURRENT_SO_VERSION} - ${LIBINT_AGE_SO_VERSION}")
message(STATUS "SO Version: Full ${LIBINT_SO_VERSION} Major ${LIBINT_MAJOR_SO_VERSION}")
