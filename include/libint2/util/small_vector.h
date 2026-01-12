/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint library.
 *
 *  Libint library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_include_libint2_util_smallvector_h_
#define _libint2_include_libint2_util_smallvector_h_

// Boost.Container's small_vector is interoperable with std::vector starting
// with version 1.61 or later see
// https://www.boost.org/doc/libs/1_61_0/doc/html/boost/container/small_vector.html#idp20337968-bb
#if defined(__has_include)
#if __has_include( \
    <boost/version.hpp>) && __has_include(<boost/container/small_vector.hpp>) && !defined(LIBINT2_DISABLE_BOOST_CONTAINER_SMALL_VECTOR)
#include <boost/version.hpp>  // read in version and do version check
#if defined(BOOST_VERSION)
#if (BOOST_VERSION / 100000 == 1) && ((BOOST_VERSION / 100 % 1000) >= 61)
#define LIBINT2_HAS_BOOST_CONTAINER_SMALL_VECTOR_H 1
#if !defined( \
    LIBINT2_SVECTOR_OPTIMIZED_RANK)  // user can override by defining
                                     // LIBINT2_SVECTOR_OPTIMIZED_RANK
#define LIBINT2_SVECTOR_OPTIMIZED_RANK 6
#endif
#include <boost/container/small_vector.hpp>
#endif  // boost version >= 1.61
#endif  // defined(BOOST_VERSION)
#else
#include <vector>
#endif
#else  // defined(__has_include)
#include <vector>
#endif  // defined(__has_include)

namespace libint2 {

#if defined(LIBINT2_HAS_BOOST_CONTAINER_SMALL_VECTOR_H)
#define LIBINT2_USES_BOOST_CONTAINER_SMALL_VECTOR_AS_SVECTOR 1
template <typename T>
using svector =
    boost::container::small_vector<T, LIBINT2_SVECTOR_OPTIMIZED_RANK>;
#else
template <typename T>
using svector = std::vector<T>;
#endif

}  // namespace libint2

#endif  // header guard
