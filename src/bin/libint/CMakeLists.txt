cmake_minimum_required(VERSION 3.8)
cmake_policy(SET CMP0074 NEW)
project(generator LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${Libint2_ROOT}/export/cmake/modules)

#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Boost 1.57 REQUIRED)

# <<<  Build  >>>

if(MSVC)
    # MSVC does not include <cmath> constants, unless _USE_MATH_DEFINES is defined.
    add_definitions("/D_USE_MATH_DEFINES")
    # Set the exception handling model
    add_definitions("/EHsc")
endif()

add_library(INT
    algebra.cc
    buildtest.cc
    class_registry.cc
    code.cc
    codeblock.cc
    comp_deriv_gauss.cc
    comp_xyz.cc
    context.cc
    default_params.cc
    dg.cc
    dgarc.cc
    dgvertex.cc
    dims.cc
    drtree.cc
    extract.cc
    flop.cc
    gauss.cc
    graph_registry.cc
    iface.cc
    iter.cc
    memory.cc
    multipole.cc
    oper.cc
    policy.cc
    policy_spec.cc
    prefactors.cc
    purgeable.cc
    rr.cc
    strategy.cc
    tactic.cc
    task.cc
    util.cc
)

add_executable(build_libint build_libint.cc)
target_link_libraries(build_libint INT Boost::boost)

add_executable(test test.cc)
target_link_libraries(test INT)

target_compile_definitions(INT PRIVATE -D__COMPILING_LIBINT2=1)
target_compile_features(INT PUBLIC "cxx_std_11")

target_include_directories(INT PUBLIC ${Libint2_BINARY_DIR}/include/
                                      ${PROJECT_SOURCE_DIR}
                                      ${Libint2_ROOT}/include/)

set_target_properties(build_libint test
    PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

#if (ENABLE_MPFR)
find_package(MPFR REQUIRED)
#get_property(_loc TARGET MPFR::MPFRXX PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
#message("iid ${_loc}")
#get_property(_loc TARGET MPFR::MPFRXX PROPERTY INTERFACE_LINK_LIBRARIES)
#message("ill ${_loc}")
#if (TARGET MPFR::GMPXX)
message("gen linking")
target_link_libraries(INT PUBLIC MPFR::MPFRXX)
message("MPFR   ${MPFR_INCLUDE} ${MPFR_LIBRARY}")
message("GMP    ${GMP_INCLUDE} ${GMP_LIBRARY}")
message("GMPXX  ${GMPXX_INCLUDE} ${GMPXX_LIBRARY}")
#endif()

#if ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
#  set_target_properties(libint_compiler PROPERTIES COMPILE_FLAGS "-c99")
#else()
#  set_target_properties(libint_compiler PROPERTIES COMPILE_FLAGS "-std=c99")
#endif()

if(MSVC)
    # Increase stack size from 1 MB to 4 MB
    set_target_properties(build_libint PROPERTIES LINK_FLAGS "/STACK:4194304")
endif()


# <<<  Install  >>>

install(TARGETS build_libint test
        RUNTIME DESTINATION bin)
