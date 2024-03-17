cmake_policy(PUSH)
cmake_policy(SET CMP0075 NEW)  # support CMAKE_REQUIRED_LIBRARIES

include(CMakePushCheckState)

cmake_push_check_state()
# needed for #include <libint2/util/vector.h>. works for top-level but not for EP; check again after more installs
list(APPEND CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}/include;${PROJECT_SOURCE_DIR}/src/lib/libint")

if(NOT(LIBINT2_REALTYPE STREQUAL "double"))
    set(LIBINT_USER_DEFINED_REAL "${LIBINT2_REALTYPE}")
    if(NOT(LIBINT_USER_DEFINED_REAL_INCLUDES))
        set(LIBINT_USER_DEFINED_REAL_INCLUDES "")
    endif()

    check_cxx_source_compiles("
//#include <libint2/util/vector.h>
${LIBINT_USER_DEFINED_REAL_INCLUDES}
  
int main(void) {
    ${LIBINT_USER_DEFINED_REAL} x1;
    ${LIBINT_USER_DEFINED_REAL} x2 = 2.0;
    ${LIBINT_USER_DEFINED_REAL} x3 = x2;
    ${LIBINT_USER_DEFINED_REAL} x4 = x2 + x3;
    ${LIBINT_USER_DEFINED_REAL} x5 = x2 - x3;
    ${LIBINT_USER_DEFINED_REAL} x6 = x2 * x3;
    ${LIBINT_USER_DEFINED_REAL} x7 = 2 * x2;
    ${LIBINT_USER_DEFINED_REAL} x8 = x2 * 3;
    x6 += x2 * x3;
    x7 -= 3 * x3;
}
"
        _user_defined_real_compiles)

    if (NOT _user_defined_real_compiles)
        message(FATAL_ERROR "LIBINT2_REALTYPE ${LIBINT_USER_DEFINED_REAL} is not usable, perhaps extra -I directories or extra #include's are needed?")
    endif()
endif()

cmake_pop_check_state()
cmake_policy(POP)
