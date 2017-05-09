cmake_minimum_required(VERSION 3.0)
project(Libintint LANGUAGES C CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(NS Libint)
set(PN ${PROJECT_NAME})

# <<<  Build  >>>

file(GLOB sources_list "${CMAKE_BINARY_DIR}/*.cc" vrr_build.c)

    # STATIC/SHARED on below governed by BUILD_SHARED_LIBS
add_library(int ${sources_list})
target_include_directories(int PRIVATE ${CMAKE_BINARY_DIR})

set_target_properties(int PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_FPIC})
if(${BUILD_SHARED_LIBS})
    target_link_libraries(int PRIVATE ${LIBC_INTERJECT})
endif()

# <<<  Install  >>>

install(FILES ${CMAKE_BINARY_DIR}/libint.h
              ${CMAKE_BINARY_DIR}/hrr_header.h
              ${CMAKE_BINARY_DIR}/vrr_header.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libint)
install(TARGETS int
        EXPORT "${PN}Targets"
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# <<<  Export Interface  >>>

math(EXPR _max_am "${MAX_AM_ERI}+1")
target_compile_definitions(int INTERFACE USING_${NS}
                                         MAX_AM_ERI=${_max_am})
target_include_directories(int INTERFACE
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
