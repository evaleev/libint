/*
 *  Copyright (C) 2004-2021 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_bin_libint_smartptr_h_
#define _libint2_src_bin_libint_smartptr_h_

#include <type_traits>
#include <libint2/config.h>

#if HAVE_SHARED_PTR_IN_BOOST
  #include <boost/shared_ptr.hpp>
  #include <boost/enable_shared_from_this.hpp>

  // For now I'll do a cheat since templated typedefs are not standard
  // Should probably at least derive SafePtr from shared_ptr
  #define SafePtr boost::shared_ptr
  #define EnableSafePtrFromThis boost::enable_shared_from_this
  #define SafePtr_from_this shared_from_this
  using boost::const_pointer_cast;
  using boost::dynamic_pointer_cast;
  using boost::static_pointer_cast;
#else
  #include <memory>
  // For now I'll do a cheat since templated typedefs are not standard
  // Should probably at least derive SafePtr from shared_ptr
  #define SafePtr std::shared_ptr
  #define EnableSafePtrFromThis std::enable_shared_from_this
  #define SafePtr_from_this shared_from_this
  using std::const_pointer_cast;
  using std::dynamic_pointer_cast;
  using std::static_pointer_cast;
#endif

namespace libint2 {
namespace detail {
/** Can be used to determine whether a type is a SafePtr */
template <typename>
struct IsSafePtrHelper : std::false_type {};
template <typename T>
struct IsSafePtrHelper<SafePtr<T>> : std::true_type {};

template <typename T>
struct IsSafePtr : IsSafePtrHelper<typename std::remove_const<
                         typename std::remove_reference<T>::type>::type> {};
}  // namespace detail
}  // namespace libint2

#endif

