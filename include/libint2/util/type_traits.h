/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2015 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License, version 2,
 *  as published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

#ifndef _libint2_include_libint2_typetraits_h_
#define _libint2_include_libint2_typetraits_h_

namespace libint2 {

  template <typename T>
  struct is_vector {
      static const bool value = false;
  };

  template <typename T>
  struct vector_traits {
      typedef T value_type;
      static const size_t extent = 1;
  };

} // namespace libint2

#endif /* _libint2_include_libint2_typetraits_h_ */
