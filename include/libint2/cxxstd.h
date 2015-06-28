/*                                                                                                        *  This file is a part of Libint.                                                                        *  Copyright (C) 2004-2014 Edward F. Valeev                                                              *                                                                                                        *  This program is free software: you can redistribute it and/or modify                                  *  it under the terms of the GNU Library General Public License, version 2,                              *  as published by the Free Software Foundation.                                                         *                                                                                                        *  This program is distributed in the hope that it will be useful,                                       *  but WITHOUT ANY WARRANTY; without even the implied warranty of                                        *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                         *  GNU General Public License for more details.                                                          *                                                                                                        *  You should have received a copy of the GNU Library General Public License                             *  along with this program.  If not, see http://www.gnu.org/licenses/.                                   *                                                                                                        */

#ifndef _libint2_include_libint2_cxxstd_h_
#define _libint2_include_libint2_cxxstd_h_

#ifndef __cplusplus
# error "Libint2 requires a C++ compiler"
#endif

#if __cplusplus >= 201402L
# define LIBINT2_CPLUSPLUS_STD 2014
#elif __cplusplus >= 201103L
# define LIBINT2_CPLUSPLUS_STD 2011
#elif __cplusplus >= 199711L
# define LIBINT2_CPLUSPLUS_STD 1998
#else
# define LIBINT2_CPLUSPLUS_STD 0 // unknown standard
#endif

// workaround: standard Intel compiler (not INDE) is not standard conforming
#if defined(__INTEL_COMPILER) && LIBINT2_CPLUSPLUS_STD==0
# ifdef __INTEL_CXX11_MODE__
#  undef LIBINT2_CPLUSPLUS_STD
#  define LIBINT2_CPLUSPLUS_STD 2011
# endif
#endif

#endif /* header guard */
