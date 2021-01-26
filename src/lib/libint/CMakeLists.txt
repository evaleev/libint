cmake_minimum_required(VERSION 3.8)
cmake_policy(SET CMP0074 NEW)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_BINARY_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/export/cmake/modules)

include(GNUInstallDirs)
include(int_computed)  # for macros.tex

# <<<  Generate Library  >>>

file(MAKE_DIRECTORY "${EXPORT_STAGE_DIR}/src")
add_custom_target(libint-generate
        COMMAND "${PROJECT_BINARY_DIR}/src/bin/libint/build_libint"
        WORKING_DIRECTORY "${EXPORT_STAGE_DIR}/src"
        COMMENT "Generating Libint2 library")
add_dependencies(libint-generate build_libint)

# <<<  Add Metadata To The Library >>>

configure_file("${PROJECT_SOURCE_DIR}/doc/progman/macros.tex.in" "${EXPORT_STAGE_DIR}/doc/macros.tex" @ONLY)
set(DATADIR_ABSOLUTE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/libint/${LIBINT_VERSION})
configure_file("${PROJECT_SOURCE_DIR}/include/libint2/basis.h.in" "${EXPORT_STAGE_DIR}/include/libint2/basis.h" @ONLY)

add_custom_target(libint-populate)
add_custom_command(TARGET libint-populate PRE_BUILD
        COMMAND ${CMAKE_COMMAND}
        "-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
        "-DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}"
        "-DEXPORT_STAGE_DIR=${EXPORT_STAGE_DIR}"
        "-DLIBINT_VERSION=${LIBINT_VERSION}"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/populate.cmake)

# <<<  Export The Library >>>

add_custom_target(libint-export
        COMMAND ${CMAKE_COMMAND} -E tar "cfvz" "${EXPORT_STAGE_DIR}.tgz" "${EXPORT_STAGE_DIR}"
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Packaging Libint2 library")
add_dependencies(libint-export libint-generate libint-populate)

