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

#ifndef _libint2_include_libint2_util_cxxstd_h_
#define _libint2_include_libint2_util_cxxstd_h_

#ifndef __cplusplus
#error "Libint2 requires a C++ compiler"
#endif

#if __cplusplus >= 202302L
#define LIBINT2_CPLUSPLUS_STD 2023
#elif __cplusplus >= 202002L
#define LIBINT2_CPLUSPLUS_STD 2020
#elif __cplusplus >= 201703L
#define LIBINT2_CPLUSPLUS_STD 2017
#elif __cplusplus >= 201402L
#define LIBINT2_CPLUSPLUS_STD 2014
#elif __cplusplus >= 201103L
#define LIBINT2_CPLUSPLUS_STD 2011
#elif __cplusplus >= 199711L
#define LIBINT2_CPLUSPLUS_STD 1998
#else
#define LIBINT2_CPLUSPLUS_STD 0  // unknown standard
#endif

// workaround: standard Intel compiler (not INDE) is not standard conforming
#if defined(__INTEL_COMPILER) && LIBINT2_CPLUSPLUS_STD == 0
#ifdef __INTEL_CXX11_MODE__
#undef LIBINT2_CPLUSPLUS_STD
#define LIBINT2_CPLUSPLUS_STD 2011
#endif
#endif

#endif /* header guard */
