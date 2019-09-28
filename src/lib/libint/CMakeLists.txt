cmake_minimum_required(VERSION 3.0)
project(library LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(NS Libint2)
set(PN ${PROJECT_NAME})

set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "ON")

# <<<  Build  >>>

#target_include_directories (chemps2-base PRIVATE ${CheMPS2_SOURCE_DIR}/CheMPS2/include/chemps2
#                                                 $<TARGET_PROPERTY:tgt::hdf5,INTERFACE_INCLUDE_DIRECTORIES>)

file(GLOB sources_list "${CMAKE_BINARY_DIR}/*.cc")

add_library(int-base OBJECT ${sources_list})
#target_include_directories(int PRIVATE ${CMAKE_BINARY_DIR})

target_compile_definitions(int-base PRIVATE -D__COMPILING_LIBINT2=1)

#target_include_directories(int-base PUBLIC ${CMAKE_BINARY_DIR}
target_include_directories(int-base PRIVATE ${CMAKE_BINARY_DIR}
                                            ${Libint2_BINARY_DIR}/include/
                                            #${PROJECT_SOURCE_DIR}
                                            ${Libint2_ROOT}/include/
                                      
)

message("which? ${BUILD_SHARED} ${BUILD_STATIC}")
if (BUILD_SHARED)
    set_target_properties(int-base PROPERTIES POSITION_INDEPENDENT_CODE 1)
endif()

if (BUILD_SHARED)
    add_library           (int-shared SHARED  $<TARGET_OBJECTS:int-base>)
    #target_link_libraries (int-shared PRIVATE ${LIBC_INTERJECT})
    #target_link_libraries (int-shared PUBLIC  tgt::lapack
    #                                          tgt::hdf5)
    set_target_properties (int-shared PROPERTIES SOVERSION ${LIBINT_MAJOR_SOVERSION}
                                                 MACOSX_RPATH ON
                                                 OUTPUT_NAME "int2"
                                                 EXPORT_NAME "int2")
    message("define shared")
endif()

if (BUILD_STATIC)
    add_library           (int-static STATIC  $<TARGET_OBJECTS:int-base>)
    #target_link_libraries (int-static PRIVATE ${LIBC_INTERJECT})
    #target_link_libraries (int-static PUBLIC  tgt::lapack
    #                                          tgt::hdf5)
    set_target_properties (int-static PROPERTIES OUTPUT_NAME "int2"
                                                 EXPORT_NAME "int2")
    message("define static")
endif()

if (BUILD_SHARED)
    add_library (int-library ALIAS int-shared)
    message("dominant shared")
elseif (BUILD_STATIC)
    add_library (int-library ALIAS int-static)
    message("dominant static")
endif()


# <<<  Install  >>>

#install(FILES ${CMAKE_BINARY_DIR}/libint.h
#              ${CMAKE_BINARY_DIR}/hrr_header.h
#              ${CMAKE_BINARY_DIR}/vrr_header.h
#        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libint)

if (BUILD_SHARED)
    install (TARGETS int-shared
             EXPORT "${PROJECT_NAME}Targets-shared"
             LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

if (BUILD_STATIC)
    install (TARGETS int-static
             EXPORT "${PROJECT_NAME}Targets-static"
             ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

## <<<  Export Interface  >>>
#
#math(EXPR _max_am "${MAX_AM_ERI}+1")
#target_compile_definitions(int INTERFACE USING_${NS}
#                                         MAX_AM_ERI=${_max_am})
#target_include_directories(int INTERFACE
#                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
#
## <<<  Export Config  >>>
#
#    # explicit "share" not "DATADIR" for CMake search path
#set(CMAKECONFIG_INSTALL_DIR "share/cmake/${NS}")
#configure_package_config_file(${PN}Config.cmake.in
#                              "${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake"
#                              INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
#        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
#install(EXPORT "${PN}Targets"
#        NAMESPACE "${NS}::"
#        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
