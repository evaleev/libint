cmake_policy(SET CMP0048 NEW)  # project_VERSION* variables populated from project(... VERSION x.x.x) string

project(Libint
        VERSION 1.2.0
        LANGUAGES C CXX)
set(Libint_AUTHORS      "Edward F. Valeev and Justin T. Fermann")
set(Libint_DESCRIPTION  "Library for the evaluation of molecular integrals of two-body operators over Gaussian functions")
set(Libint_URL          "http://libint.valeyev.net")
set(Libint_LICENSE      "GPL-3.0 for generator; LGPL-3.0 for generated")

cmake_minimum_required(VERSION 3.0)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

################################### Options ####################################
include(psi4OptionsTools)
option_with_default(CMAKE_BUILD_TYPE "Build type" Release)
option_with_default(MAX_AM_ERI "Maximum angular momentum for integrals (min value 3)" 5)
option_with_print(BUILD_SHARED_LIBS "Build final library as shared, not static" OFF)
option_with_default(BUILD_FPIC "Libraries will be compiled with position independent code" ON)
if(${BUILD_SHARED_LIBS} AND NOT ${BUILD_FPIC})
    message(FATAL_ERROR "BUILD_SHARED_LIBS ON and BUILD_FPIC OFF are incompatible, as shared library requires position independent code")
endif()
option_with_print(ENABLE_GENERIC "Enable mostly static linking in shared library" OFF)
option_with_flags(ENABLE_XHOST "Enable processor-specific optimization" ON
                    "-xHost" "-march=native")
option_with_print(MERGE_LIBDERIV_INCLUDEDIR "Install libderiv headers to libint namespace. Psi4 wants ON" OFF)

######################### Process & Validate Options ###########################
include(autocmake_safeguards)
include(autocmake_static_library)

################################# Main Project #################################
set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/staged)

add_subdirectory(src)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(PN ${PROJECT_NAME})

# <<<  Install  >>>

install(DIRECTORY ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
                  ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
                  ${STAGED_INSTALL_PREFIX}/share
        DESTINATION ${CMAKE_INSTALL_PREFIX}
        USE_SOURCE_PERMISSIONS
        PATTERN "*header*" EXCLUDE)

# <<<  Export Config  >>>

# GNUInstallDirs "DATADIR" wrong here; CMake search path wants "share".
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PN}")
configure_package_config_file(cmake/${PN}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake"
                              INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
                                 VERSION ${${PN}_VERSION}
                                 COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
