cmake_minimum_required(VERSION 3.8)
cmake_policy(SET CMP0074 NEW)
project(generator LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${Libint2_ROOT}/export/cmake/modules)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# <<<  Build  >>>

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
target_link_libraries(build_libint INT)

add_executable(test test.cc)
target_link_libraries(test INT)

target_compile_definitions(INT
    PRIVATE -D__COMPILING_LIBINT2=1)

target_include_directories(INT PUBLIC ${PROJECT_SOURCE_DIR}
                                      ${Libint2_ROOT}/include/
                                      ${Libint2_BINARY_DIR}/include/
)

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
message("MPFR ${MPFR_INCLUDE}")
message("GMP  ${GMPXX_INCLUDE}")

#if ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
#  set_target_properties(libint_compiler PROPERTIES COMPILE_FLAGS "-c99")
#else()
#  set_target_properties(libint_compiler PROPERTIES COMPILE_FLAGS "-std=c99")
#endif()
#
#if(MSVC)
#    # Increase stack size from 1 MB to 4 MB
#    set_target_properties(libint_compiler PROPERTIES LINK_FLAGS "/STACK:4194304")
#endif()

# <<<  Install  >>>

install(TARGETS build_libint test
        RUNTIME DESTINATION bin)
