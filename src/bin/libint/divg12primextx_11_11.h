
#ifndef _libint2_src_bin_libint_divg12primextx1111_h_
#define _libint2_src_bin_libint_divg12primextx1111_h_

#include <integral.h>
#include <integral_11_11.h>

using namespace std;

namespace libint2 {
  
  /**
     DivG12prime_xTx_11_11 --
     integral over DivG12prime_xTx operator with one bfs for each particle in bra and ket.
  */
  typedef GenIntegralSet_11_11<DivG12prime_xTx,CGShell,EmptySet> DivG12prime_xTx_11_11_sq;
  typedef GenIntegralSet_11_11<DivG12prime_xTx,CGF,EmptySet> DivG12prime_xTx_11_11_int;

};

#endif

