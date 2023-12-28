add_library(
  libint-libcompiler
  STATIC
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

target_compile_definitions(
  libint-libcompiler
  PUBLIC
    __COMPILING_LIBINT2=1
    # MSVC does not include <cmath> constants, unless _USE_MATH_DEFINES is defined.
    $<$<BOOL:${MSVC}>:_USE_MATH_DEFINES>
  )

target_compile_options(
  libint-libcompiler
  PUBLIC
    # Set the exception handling model (allows "throw")
    $<$<BOOL:${MSVC}>:"/EHsc">
  )

target_compile_features(
  libint-libcompiler
  PUBLIC
    "cxx_std_11"
  )

set_target_properties(
  libint-libcompiler
  PROPERTIES
    CXX_EXTENSIONS OFF
    UNITY_BUILD FALSE  # do not use unity build for build_libint
    LIBRARY_OUTPUT_NAME int2-libcompiler
    ARCHIVE_OUTPUT_NAME int2-libcompiler
  )

target_include_directories(
  libint-libcompiler
  PUBLIC
    ${PROJECT_BINARY_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/include
  )

target_link_libraries(
  libint-libcompiler
  PUBLIC
    Boost::headers
    Multiprecision::gmpxx
  )

add_executable(
  build_libint
  EXCLUDE_FROM_ALL
  build_libint.cc
  )
target_link_libraries(build_libint libint-libcompiler)

if(MSVC)
    # Increase stack size from 1 MB to 4 MB
    set_target_properties(
      build_libint
      PROPERTIES
        LINK_OPTIONS "/STACK:4194304"
      )
endif()

set_target_properties(
  build_libint
  PROPERTIES
    INSTALL_RPATH_USE_LINK_PATH TRUE
  )
