# <<<  Build  >>>

if(MSVC)
    # MSVC does not include <cmath> constants, unless _USE_MATH_DEFINES is defined.
    add_compile_definitions(_USE_MATH_DEFINES)
    # Set the exception handling model
    add_compile_options(/EHsc)
endif()

### compiler library

add_library(libint-libcompiler STATIC
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

target_compile_definitions(libint-libcompiler PUBLIC -D__COMPILING_LIBINT2=1)
target_compile_features(libint-libcompiler PUBLIC "cxx_std_11")
set_target_properties(libint-libcompiler PROPERTIES CXX_EXTENSIONS OFF)
set_target_properties(libint-libcompiler PROPERTIES UNITY_BUILD FALSE) # do not use unity build for build_libint

target_include_directories(libint-libcompiler PUBLIC ${PROJECT_BINARY_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/include)
set_target_properties(libint-libcompiler PROPERTIES LIBRARY_OUTPUT_NAME int2-libcompiler ARCHIVE_OUTPUT_NAME int2-libcompiler)

## compiler library dependencies

find_package(Boost 1.57 REQUIRED)
find_package(Multiprecision MODULE REQUIRED COMPONENTS gmpxx)
target_link_libraries(
  libint-libcompiler
  PUBLIC
    Boost::headers
    Multiprecision::gmpxx
  )

get_property(_loc TARGET Multiprecision::gmp PROPERTY LOCATION)
message(STATUS "${Cyan}Found GMP${ColourReset}: ${_loc}")
get_property(_loc TARGET Multiprecision::gmpxx PROPERTY LOCATION)
message(STATUS "${Cyan}Found GMPXX${ColourReset}: ${_loc}")
if (TARGET Multiprecision::mpfr)
    get_property(_loc TARGET Multiprecision::mpfr PROPERTY LOCATION)
    message(STATUS "${Cyan}Found MPFR${ColourReset}: ${_loc} (found version ${MPFR_VERSION})")
endif()

### executables compiler and test

add_executable(build_libint EXCLUDE_FROM_ALL build_libint.cc)
target_link_libraries(build_libint libint-libcompiler)
if(MSVC)
    # Increase stack size from 1 MB to 8 MB
    set_target_properties(build_libint PROPERTIES LINK_FLAGS "/STACK:8388608")
endif()

set_target_properties(build_libint
    PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

