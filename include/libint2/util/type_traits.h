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

#ifndef _libint2_include_libint2_typetraits_h_
#define _libint2_include_libint2_typetraits_h_

#include <cstdlib>

namespace libint2 {

template <typename T>
struct is_vector {
  static const bool value = false;
};

template <typename T>
struct vector_traits {
  typedef T scalar_type;
  static const std::size_t extent = 1;
};

}  // namespace libint2

#endif /* _libint2_include_libint2_typetraits_h_ */
