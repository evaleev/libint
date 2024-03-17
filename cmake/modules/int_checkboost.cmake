cmake_policy(PUSH)
cmake_policy(SET CMP0075 NEW)  # support CMAKE_REQUIRED_LIBRARIES

include(CMakePushCheckState)

cmake_push_check_state()
list(APPEND CMAKE_REQUIRED_LIBRARIES Boost::headers)
#list(APPEND CMAKE_REQUIRED_FLAGS "-std=c++11")  # set CMAKE_CXX_STANDARD and 0067 NEW ?

check_cxx_source_compiles("
#include <boost/preprocessor.hpp>

int main(void) {  
#if not BOOST_PP_VARIADICS  // no variadic macros? your compiler is out of date! (should not be possible since variadic macros are part of C++11)
#  error \"your compiler does not provide variadic macros (but does support C++11), something is seriously broken, please create an issue at https://github.com/evaleev/libint/issues\"
#endif
    return 0;
}
"
    _boost_pp_variadics)

if (NOT _boost_pp_variadics)
    message(FATAL_ERROR "BOOST_PP_VARIADICS is oddly missing from detected installation")
endif()

cmake_pop_check_state()
cmake_policy(POP)
