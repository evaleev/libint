/*
  Define intrinsic data types here
*/

#ifndef _libint2_include_libint2intrinsictypes_h_
#define _libint2_include_libint2intrinsictypes_h_

#include <libint2_config.h>

/* default LIBINT2 real type is double, but user can specify a different type to use with configure */
#ifdef LIBINT2_USER_DEFINED_REAL
typedef LIBINT2_USER_DEFINED_REAL LIBINT2_FLOAT;
#else
typedef double LIBINT2_FLOAT;
#endif

#include <limits.h>

/* determine default LIBINT2 64-bit integer */
#ifdef HAVE_STDINT_H

#include <stdint.h>
/* because mpz_class does not mesh with long long types, only use those when absolutely necessary */
#if UINT_LEAST64_MAX != ULONG_MAX
  typedef uint_least64_t LIBINT2_UINT_LEAST64;
#else
  typedef unsigned long int LIBINT2_UINT_LEAST64;
#endif
#if INT_LEAST64_MAX != LONG_MAX
  typedef int_least64_t LIBINT2_INT_LEAST64;
#else
  typedef long int LIBINT2_INT_LEAST64;
#endif

#else

#if defined(ULONGLONG_MAX) && !defined(ULLONG_MAX)
#    define ULLONG_MAX ULONGLONG_MAX
#endif

# ifdef ULLONG_MAX
#   if ULONGLONG_MAX == (0xffffffffffffffffuLL) /* uLL reqd for xlC */
     typedef long long LIBINT2_INT_LEAST64; 
     typedef unsigned long long LIBINT2_UINT_LEAST64; 
#   else
#     error defaults not correct; you must hand modify libint2_intrinsic_types.h
#   endif
# elif ULONG_MAX != 0xffffffff

#   if ULONG_MAX == 18446744073709551615 /* 2**64 - 1 */
     typedef long LIBINT2_INT_LEAST64;
     typedef unsigned long LIBINT2_UINT_LEAST64;
#   else
#     error defaults not correct; you must hand modify scint.h
#   endif
# else /* assume no 64-bit integers */
#   error 64 bit integer types are required
# endif

#endif /* HAVE_STDINT_H */

#endif
