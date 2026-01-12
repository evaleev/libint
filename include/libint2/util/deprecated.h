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

#ifndef _libint2_include_libint2_util_deprecated_h_
#define _libint2_include_libint2_util_deprecated_h_

#ifndef LIBINT2_DEPRECATED  // avoid clashing with previous definitions
#if __cplusplus >= 201402L
#define LIBINT2_DEPRECATED [[deprecated]]
#elif defined(__GNUC__)
#define LIBINT2_DEPRECATED __attribute__((deprecated))
#else
#pragma message("WARNING: You need to implement LIBINT2_DEPRECATED for this compiler")
#define LIBINT2_DEPRECATED
#endif
#endif  // not defined(LIBINT2_DEPRECATED)

#endif /* header guard */
