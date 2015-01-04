// shamelessly borrowed from MADNESS (GPLv3)
// written by Justus Calvin (justus.c79@gmail.com)

#ifndef _libint2_include_stdarray_h
#define _libint2_include_stdarray_h

#include <libint2/config.h>

#if defined(LIBINT_USE_ARRAY)
#  include <array>
#elif defined(LIBINT_USE_TR1_ARRAY)
#  include <tr1/array>
#elif defined(LIBINT_USE_BOOST_TR1_ARRAY_HPP)
#  include <boost/tr1/array.hpp>
#else
#  define LIBINT_HAS_STD_ARRAY 1
#  include <stdarray_bits.h>
   namespace std {
      using namespace libint2::tr1::array;
   }
#endif

namespace std {
#if defined(LIBINT_HAS_STD_TR1_ARRAY) && !defined(LIBINT_HAS_STD_ARRAY)
#   define LIBINT_HAS_STD_ARRAY 1
    using ::std::tr1::array;
    using ::std::tr1::swap;
    using ::std::tr1::tuple_size;
    using ::std::tr1::tuple_element;
    using ::std::tr1::tuple_size;
    using ::std::tr1::tuple_element;
    using ::std::tr1::get;
#endif
}

#endif // header guard
