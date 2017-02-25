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

#ifndef _libint2_src_lib_libint_libint2memory_h_
#define _libint2_src_lib_libint_libint2memory_h_

#include <cstdlib>
#include <libint2/libint2_params.h>

namespace libint2 {

  /// Aligned version of malloc().

  /** Allocates a memory block aligned to LIBINT2_ALIGN_SIZE*sizeof(LIBINT2_REALTYPE).
      If LIBINT2_ALIGN_SIZE, no alignment is assumed. Use free() to deallocate. */
  inline void* malloc(size_t nbytes) {
    void* result;
#if (LIBINT2_ALIGN_SIZE == 0)
    result = ::malloc(nbytes);
#elif defined(HAVE_POSIX_MEMALIGN)
    posix_memalign(&result, LIBINT2_ALIGN_SIZE*sizeof(LIBINT2_REALTYPE), nbytes);
#else
#   error "LIBINT2_ALIGN_SIZE!=0 but posix_memalign is not available"
#endif
    return result;
  }

  /// type-specific version of libint2::malloc()
  template <typename T>
  inline T* malloc(size_t n) {
    return reinterpret_cast<T*>(malloc(n * sizeof(T)));
  }

}

#endif /* _libint2_src_lib_libint_libint2memory_h_ */
