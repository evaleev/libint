cmake_minimum_required(VERSION 3.0)
project(library LANGUAGES CXX)

#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(NS Libint2)
set(PN ${PROJECT_NAME})

find_package(Eigen3 REQUIRED)
    #get_property(_ill TARGET Eigen3::Eigen PROPERTY INTERFACE_LINK_LIBRARIES)
    #get_property(_iid TARGET Eigen3::Eigen PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
    #message("ill ${_ill}")
    #message("iid ${_iid}")


# <<<  Build  >>>

#target_include_directories (chemps2-obj PRIVATE ${CheMPS2_SOURCE_DIR}/CheMPS2/include/chemps2
#                                                 $<TARGET_PROPERTY:tgt::hdf5,INTERFACE_INCLUDE_DIRECTORIES>)

file(GLOB sources_list "${CMAKE_BINARY_DIR}/*.cc")

add_library(int-obj OBJECT ${sources_list})
#target_include_directories(int PRIVATE ${CMAKE_BINARY_DIR})

target_compile_definitions(int-obj PRIVATE -D__COMPILING_LIBINT2=1)
target_compile_features(int-obj PUBLIC "cxx_std_11")

#target_include_directories(int-obj PUBLIC ${CMAKE_BINARY_DIR}
target_include_directories(int-obj PRIVATE ${CMAKE_BINARY_DIR}
                                            ${Libint2_BINARY_DIR}/include/
                                            ${CMAKE_SOURCE_DIR}
                                            ${Libint2_ROOT}/include/
)

if (BUILD_SHARED)
    set_target_properties(int-obj PROPERTIES POSITION_INDEPENDENT_CODE 1)
endif()

if (BUILD_SHARED)
    add_library           (int-shared SHARED  $<TARGET_OBJECTS:int-obj>)
    target_link_libraries (int-shared INTERFACE Eigen3::Eigen)
    set_target_properties (int-shared PROPERTIES SOVERSION ${LIBINT_MAJOR_SOVERSION}
                                                 MACOSX_RPATH ON
                                                 OUTPUT_NAME "int2"
                                                 EXPORT_NAME "cxx")
    target_include_directories(int-shared INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
endif()

if (BUILD_STATIC)
    add_library           (int-static STATIC  $<TARGET_OBJECTS:int-obj>)
    target_link_libraries (int-static INTERFACE Eigen3::Eigen)
    set_target_properties (int-static PROPERTIES OUTPUT_NAME "int2"
                                                 EXPORT_NAME "cxx")
    target_include_directories(int-static INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
endif()

if (BUILD_SHARED)
    add_library (int-library ALIAS int-shared)
elseif (BUILD_STATIC)
    add_library (int-library ALIAS int-static)
endif()


# <<<  Install  >>>

#install(FILES ${CMAKE_BINARY_DIR}/libint.h
#              ${CMAKE_BINARY_DIR}/hrr_header.h
#              ${CMAKE_BINARY_DIR}/vrr_header.h
install(DIRECTORY ${Libint2_BINARY_DIR}/include/
                  ${Libint2_ROOT}/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PATTERN "*.in$" EXCLUDE)
        #PATTERN "*.h$"
        #PATTERN "*.hpp$")
install(FILES ${CMAKE_BINARY_DIR}/libint2_iface.h
              ${CMAKE_BINARY_DIR}/libint2_iface_internal.h
              ${CMAKE_BINARY_DIR}/libint2_params.h
              ${CMAKE_BINARY_DIR}/libint2_types.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libint2)

target_include_directories(int-obj PRIVATE ${CMAKE_BINARY_DIR})

if (BUILD_SHARED)
    install (TARGETS int-shared
             EXPORT "${NS}Targets-shared"
             LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

if (BUILD_STATIC)
    install (TARGETS int-static
             EXPORT "${NS}Targets-static"
             ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# <<<  Export Interface  >>>

#math(EXPR _max_am "${MAX_AM_ERI}+1")
#target_compile_definitions(int INTERFACE USING_${NS}
#                                         MAX_AM_ERI=${_max_am})

# <<<  Export Config  >>>

# explicit "share" not "DATADIR" for CMake search path
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${NS}")
#configure_package_config_file(${PN}Config.cmake.in
#                              "${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake"
#                              INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
#        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

if (BUILD_SHARED)
    install(EXPORT "${NS}Targets-shared"
                   NAMESPACE "${NS}::"
                   DESTINATION ${CMAKECONFIG_INSTALL_DIR})
endif()

if (BUILD_STATIC)
    install(EXPORT "${NS}Targets-static"
                   NAMESPACE "${NS}::"
                   DESTINATION ${CMAKECONFIG_INSTALL_DIR})
endif()
