cmake_minimum_required(VERSION 3.0)
project(libint_compiler LANGUAGES C)

# <<<  Build  >>>

# create dummy header file
configure_file(../libint_config.h.cmake.in libint_config.h @ONLY)

set(COMPILER_SRC build_libint.c
                 emit_hrr_build.c
                 emit_hrr_build_macro.c
                 emit_order.c
                 emit_vrr_build.c
                 emit_vrr_build_macro.c
                 mem_man.c
)
add_executable(libint_compiler ${COMPILER_SRC})

# link -lm only if necessary
find_package(StandardMathLibraryC)
target_link_libraries(libint_compiler ${STANDARD_MATH_LIBRARY})
target_include_directories(libint_compiler PRIVATE ../
                                                   ${PROJECT_BINARY_DIR})

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
  set_target_properties(libint_compiler PROPERTIES COMPILE_FLAGS "-c99")
else()
  set_target_properties(libint_compiler PROPERTIES COMPILE_FLAGS "-std=c99")
endif()

if(MSVC)
    # Increase stack size from 1 MB to 4 MB
    set_target_properties(libint_compiler PROPERTIES LINK_FLAGS "/STACK:4194304")
endif(MSVC)

math(EXPR LIBINT_MAX_AM ${MAX_AM_ERI}+1)
math(EXPR LIBINT_OPT_AM "${MAX_AM_ERI}/2 + ${MAX_AM_ERI}%2")
target_compile_definitions(libint_compiler
    PRIVATE -DLIBINT_MAX_AM=${LIBINT_MAX_AM}
            -DLIBINT_OPT_AM=${LIBINT_OPT_AM}
            -DLIBINT_MAX_CLASS_SIZE=300
            -DLONG_DOUBLE=0)

# <<<  Install  >>>

install(TARGETS libint_compiler
        RUNTIME DESTINATION bin)
