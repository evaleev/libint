cmake_minimum_required(VERSION 3.0)
project(libderiv_compiler LANGUAGES C)

# <<<  Build  >>>

# create dummy header file
configure_file(../libint_config.h.cmake.in libint_config.h @ONLY)

set(COMPILER_SRC build_libderiv.c
                 emit_d1hrr_build_macro.c
                 emit_deriv1_managers.c
                 emit_deriv_build_macro.c
                 emit_d1hrr_build.c
                 emit_deriv12_managers.c
                 emit_deriv_build.c
                 mem_man.c
)
add_executable(libderiv_compiler ${COMPILER_SRC})

target_include_directories(libderiv_compiler PRIVATE ../
                                                     ${PROJECT_BINARY_DIR})

math(EXPR LIBDERIV_OPT_AM1 ${MAX_AM_ERI}-1)  # A.M. level for 1st derivative ERIs
math(EXPR LIBINT_MAX_AM ${MAX_AM_ERI}+1)
math(EXPR LIBINT_OPT_AM "${MAX_AM_ERI}/2 + ${MAX_AM_ERI}%2")
target_compile_definitions(libderiv_compiler
    PRIVATE -DLIBDERIV_MAX_AM1=${MAX_AM_ERI}
            -DLIBDERIV_MAX_AM12=${LIBDERIV_OPT_AM1}
            -DLIBINT_MAX_AM=${LIBINT_MAX_AM}
            -DLIBINT_OPT_AM=${LIBINT_OPT_AM})

# acquire libint.h
find_package(Libintint CONFIG REQUIRED)
target_include_directories(libderiv_compiler PRIVATE
    $<TARGET_PROPERTY:Libint::int,INTERFACE_INCLUDE_DIRECTORIES>)

# <<<  Install  >>>

install(TARGETS libderiv_compiler
        RUNTIME DESTINATION bin)
