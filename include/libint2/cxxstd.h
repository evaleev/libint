/*                                                                                                        *  This file is a part of Libint.                                                                        *  Copyright (C) 2004-2014 Edward F. Valeev                                                              *                                                                                                        *  This program is free software: you can redistribute it and/or modify                                  *  it under the terms of the GNU Library General Public License, version 2,                              *  as published by the Free Software Foundation.                                                         *                                                                                                        *  This program is distributed in the hope that it will be useful,                                       *  but WITHOUT ANY WARRANTY; without even the implied warranty of                                        *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                         *  GNU General Public License for more details.                                                          *                                                                                                        *  You should have received a copy of the GNU Library General Public License                             *  along with this program.  If not, see http://www.gnu.org/licenses/.                                   *                                                                                                        */

#ifndef _libint2_include_libint2_cxxstd_h_
#define _libint2_include_libint2_cxxstd_h_

#ifndef __cplusplus
# error "Libint2 requires a C++ compiler"
#endif

#if __cplusplus >= 201103L
#  define _libint2_has_cpp11_
#endif

// workaround: Intel compiler is not standard conforming

#ifdef __INTEL_COMPILER
# ifdef __INTEL_CXX11_MODE__
#   define _libint2_has_cpp11_
# endif
#endif

#if not _libint2_has_cpp11_
# error " Libint2 C++ API requires C++11 support"
#endif

#endif /* header guard */

