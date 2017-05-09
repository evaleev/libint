cmake_minimum_required(VERSION 3.0)
project(Libintderiv LANGUAGES C CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(NS Libint)
set(PN ${PROJECT_NAME})

# <<<  Build  >>>

file(GLOB sources_list "${CMAKE_BINARY_DIR}/*.cc" vrr_build.c)

    # STATIC/SHARED on below governed by BUILD_SHARED_LIBS
add_library(deriv ${sources_list})
target_include_directories(deriv PRIVATE ${CMAKE_BINARY_DIR})
    # acquire libint.h and libint library
find_package(Libintint CONFIG REQUIRED)
target_link_libraries(deriv PRIVATE Libint::int)

set_target_properties(deriv PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_FPIC})
if(${BUILD_SHARED_LIBS})
    target_link_libraries(deriv PRIVATE ${LIBC_INTERJECT})
endif()

# <<<  Install  >>>

if(${MERGE_LIBDERIV_INCLUDEDIR})
    install(FILES ${CMAKE_BINARY_DIR}/libderiv.h
                  ${CMAKE_BINARY_DIR}/d1hrr_header.h
                  ${CMAKE_BINARY_DIR}/deriv_header.h
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libint)
else()
    install(FILES ${CMAKE_BINARY_DIR}/libderiv.h
                  ${CMAKE_BINARY_DIR}/d1hrr_header.h
                  ${CMAKE_BINARY_DIR}/deriv_header.h
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libderiv)
endif()

install(TARGETS deriv
        EXPORT "${PN}Targets"
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# <<<  Export Interface  >>>

math(EXPR _max_am "${MAX_AM_ERI}+1")
target_compile_definitions(deriv INTERFACE USING_${NS}
                                           MAX_AM_ERI=${_max_am})
target_include_directories(deriv INTERFACE
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# <<<  Export Config  >>>

    # explicit "share" not "DATADIR" for CMake search path
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${NS}")
configure_package_config_file(${PN}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake"
                              INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(EXPORT "${PN}Targets"
        NAMESPACE "${NS}::"
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
