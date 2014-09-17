/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
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

#ifndef _libint2_src_lib_libint_cxxapi_h_
#define _libint2_src_lib_libint_cxxapi_h_

#if __cplusplus <= 199711L
# error " Libint2 C++ API requires C++11 support"
#endif

#include <libint2.h>

namespace libint2 {
  inline void init() {
    libint2_static_init();
  }
  inline void cleanup() {
    libint2_static_cleanup();
  }
}

#include <libint2/engine.h>

#endif /* _libint2_src_lib_libint_cxxapi_h_ */
