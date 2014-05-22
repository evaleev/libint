/*
 * memory.h
 *
 *  Created on: May 22, 2014
 *      Author: evaleev
 */

#ifndef _libint2_src_lib_libint_libint2memory_h_
#define _libint2_src_lib_libint_libint2memory_h_

#include <cstdlib>
#include <libint2_params.h>

namespace libint2 {

  /// Aligned version of malloc().

  /** Allocates a memory block aligned to LIBINT2_ALIGN_SIZE*sizeof(LIBINT2_REALTYPE).
      If LIBINT2_ALIGN_SIZE, no alignment is assumed Use free() to deallocate. */
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
