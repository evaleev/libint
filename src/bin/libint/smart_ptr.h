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

#include <libint2/config.h>
#include <memory>

namespace libint2 {
namespace detail {
/** Can be used to determine whether a type is a std::shared_ptr */
template <typename T>
struct IsSafePtr {
  enum { result = false };
};

template <typename T>
struct IsSafePtr< std::shared_ptr<T> > {
  enum { result = true };
};
template <typename T>
struct IsSafePtr< const std::shared_ptr<T> > {
  enum { result = true };
};
template <typename T>
struct IsSafePtr< std::shared_ptr<T>& > {
  enum { result = true };
};
template <typename T>
struct IsSafePtr< const std::shared_ptr<T>& > {
  enum { result = true };
};
}  // namespace detail
}  // namespace libint2

#endif

