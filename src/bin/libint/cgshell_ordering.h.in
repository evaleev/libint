
#ifndef _libint2_src_bin_libint_cgshellordering_h_
#define _libint2_src_bin_libint_cgshellordering_h_

#include <libint2_config.h>
/* if this is not defined, then using exported source -- redefine using macro from libint2_params.h */
#ifndef LIBINT_CARTGAUSS_MAX_AM
# include <libint2_params.h>
# define LIBINT_CARTGAUSS_MAX_AM LIBINT2_CARTGAUSS_MAX_AM
#endif
#ifndef LIBINT_CGSHELL_ORDERING
# include <libint2_params.h>
# define LIBINT_CGSHELL_ORDERING LIBINT2_CGSHELL_ORDERING
# define LIBINT_CGSHELL_ORDERING_STANDARD LIBINT2_CGSHELL_ORDERING_STANDARD
# define LIBINT_CGSHELL_ORDERING_INTV3 LIBINT2_CGSHELL_ORDERING_INTV3
# define LIBINT_CGSHELL_ORDERING_GAMESS LIBINT2_CGSHELL_ORDERING_GAMESS
# define LIBINT_CGSHELL_ORDERING_ORCA LIBINT2_CGSHELL_ORDERING_ORCA
# define LIBINT_CGSHELL_ORDERING_BAGEL LIBINT2_CGSHELL_ORDERING_BAGEL
#endif

namespace libint2 {

  enum CGShellOrdering {
    CGShellOrdering_Standard     = LIBINT_CGSHELL_ORDERING_STANDARD,
    CGShellOrdering_IntV3        = LIBINT_CGSHELL_ORDERING_INTV3,
    CGShellOrdering_GAMESS       = LIBINT_CGSHELL_ORDERING_GAMESS,
    CGShellOrdering_ORCA         = LIBINT_CGSHELL_ORDERING_ORCA,
    CGShellOrdering_BAGEL        = LIBINT_CGSHELL_ORDERING_BAGEL
  };

};

#include <cgshellinfo.h>

//
// Macros common to all orderings
//

/* Computes the number of Cartesian function in a shell given
 * am = total angular momentum
 * formula: (am*(am+1))/2 + am+1;
 */
#define INT_NCART(am) ((((am)+2)*((am)+1))>>1)

/* For a given ang. mom., am, with n cartesian functions, compute the
 * number of cartesian functions for am+1 or am-1
 */
#define INT_NCART_DEC(am,n) ((n)-(am)-1)
#define INT_NCART_INC(am,n) ((n)+(am)+2)

//
// Macros that define orderings
//

#if LIBINT_CGSHELL_ORDERING == LIBINT_CGSHELL_ORDERING_STANDARD
// this piece of code is from MPQC:src/lib/chemistry/cca/int/macros.h
// Copyright Edward Valeev

/* Computes an index to a Cartesian function within a shell given
 * am = total angular momentum
 * i = the exponent of x (i is used twice in the macro--beware side effects)
 * j = the exponent of y
 * formula: (am - i + 1)*(am - i)/2 + am - i - j unless i==am, then 0
 * The following loop will generate indices in the proper order:
 *  cartindex = 0;
 *  for (i=am; i>=0; i--) {
 *    for (j=am-i; j>=0; j--) {
 *      do_it_with(cartindex);
 *      cartindex++;
 *      }
 *    }
 */
#define INT_CARTINDEX(am,i,j) ( (((((am) - (i) + 1)*((am) - (i)))>>1) + (am) - (i) - (j)) )

/* This sets up the above loop over cartesian exponents as follows
 * int i, j, k;
 * FOR_CART(i,j,k,am)
 *   Stuff using i,j,k.
 *   END_FOR_CART
 */
#define FOR_CART(i,j,k,am) for((i)=(am);(i)>=0;(i)--) {\
                           for((j)=(am)-(i);(j)>=0;(j)--) \
                           { (k) = (am) - (i) - (j);
#define END_FOR_CART }}

#endif // STANDARD ordering

#if LIBINT_CGSHELL_ORDERING == LIBINT_CGSHELL_ORDERING_INTV3
// this piece of code is from MPQC:src/lib/chemistry/qc/intv3/macros.h
// Copyright Curtis Janssen

/* Computes an index to a Cartesian function within a shell given
 * am = total angular momentum
 * i = the exponent of x (i is used twice in the macro--beware side effects)
 * j = the exponent of y
 * formula: am*(i+1) - (i*(i+1))/2 + i+1 - j - 1
 * The following loop will generate indices in the proper order:
 *  cartindex = 0;
 *  for (i=0; i<=am; i++) {
 *    for (k=0; k<=am-i; k++) {
 *      j = am - i - k;
 *      do_it_with(cartindex); // cartindex == INT_CARTINDEX(am,i,j)
 *      cartindex++;
 *      }
 *    }
 */
#define INT_CARTINDEX(am,i,j) (((((((am)+1)<<1)-(i))*((i)+1))>>1)-(j)-1)

/* This sets up the above loop over cartesian exponents as follows
 * FOR_CART(i,j,k,am)
 *   Stuff using i,j,k.
 *   END_FOR_CART
 */
#define FOR_CART(i,j,k,am) for((i)=0;(i)<=(am);(i)++) {\
                           for((k)=0;(k)<=(am)-(i);(k)++) \
                           { (j) = (am) - (i) - (k);
#define END_FOR_CART }}

#endif // INTV3 ordering

#if LIBINT_CGSHELL_ORDERING == LIBINT_CGSHELL_ORDERING_GAMESS
// for definition of the ordering see CGShellInfo

/* Computes an index to a Cartesian function within a shell given
 * am = total angular momentum
 * i = the exponent of x (i is used twice in the macro--beware side effects)
 * j = the exponent of y
 * for this ordering there is no formula
 */
#define INT_CARTINDEX(am,i,j) CGShellInfo< CGShellOrderingData<CGShellOrdering_GAMESS,LIBINT_CARTGAUSS_MAX_AM> >::cartindex(am,i,j)

/* This sets up the above loop over cartesian exponents as follows
 * FOR_CART(i,j,k,am)
 *   Stuff using i,j,k.
 *   END_FOR_CART
 */
#define FOR_CART(i,j,k,am) for(int __xyz=0; __xyz<INT_NCART(am); ++__xyz) { \
                             CGShellInfo< CGShellOrderingData<CGShellOrdering_GAMESS,LIBINT_CARTGAUSS_MAX_AM> >::cartindex_to_ijk(am,__xyz,i,j,k);
#define END_FOR_CART }

#endif // GAMESS ordering

#if LIBINT_CGSHELL_ORDERING == LIBINT_CGSHELL_ORDERING_ORCA
// for definition of the ordering see CGShellInfo

/* Computes an index to a Cartesian function within a shell given
 * am = total angular momentum
 * i = the exponent of x (i is used twice in the macro--beware side effects)
 * j = the exponent of y
 * for this ordering there is no formula
 */
#define INT_CARTINDEX(am,i,j) CGShellInfo< CGShellOrderingData<CGShellOrdering_ORCA,LIBINT_CARTGAUSS_MAX_AM> >::cartindex(am,i,j)

/* This sets up the above loop over cartesian exponents as follows
 * FOR_CART(i,j,k,am)
 *   Stuff using i,j,k.
 *   END_FOR_CART
 */
#define FOR_CART(i,j,k,am) for(int __xyz=0; __xyz<INT_NCART(am); ++__xyz) { \
                             CGShellInfo< CGShellOrderingData<CGShellOrdering_ORCA,LIBINT_CARTGAUSS_MAX_AM> >::cartindex_to_ijk(am,__xyz,i,j,k);
#define END_FOR_CART }

#endif // ORCA ordering

#if LIBINT_CGSHELL_ORDERING == LIBINT_CGSHELL_ORDERING_BAGEL
// permuted version of IntV3 (y in IntV3 = x in Bagel, etc.)

/* Computes an index to a Cartesian function within a shell given
 * all arguments are used multiple times in the macro--beware side effects)
 * am = total angular momentum
 * i = the exponent of x
 * j = the exponent of y
 * formula: am*(k+1) - (k*(k+1))/2 + k+1 - i - 1 = 
 * The following loop will generate indices in the proper order:
 *  cartindex = 0;
 *  for (k=0; k<=am; k++) {
 *    for (j=0; j<=am-k; j++) {
 *      i = am - j - k;
 *      do_it_with(cartindex); // cartindex == INT_CARTINDEX(am,i,j)
 *      cartindex++;
 *      }
 *    }
 */
#define INT_CARTINDEX(am,i,j) ((((am)+((i)+(j))+2)*((am)-((i)+(j))+1)>>1)-(i)-1)

/* This sets up the above loop over cartesian exponents as follows
 * FOR_CART(i,j,k,am)
 *   Stuff using i,j,k.
 *   END_FOR_CART
 */
#define FOR_CART(i,j,k,am) for((k)=0;(k)<=(am);(k)++) {\
                           for((j)=0;(j)<=(am)-(k);(j)++) \
                           { (i) = (am) - (j) - (k);
#define END_FOR_CART }}

#endif // Bagel ordering

#endif // header guard
