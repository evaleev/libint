
#include <integral.h>

#ifndef _libint2_src_bin_libint_dummyintegral_h_
#define _libint2_src_bin_libint_dummyintegral_h_

namespace libint2 {
  typedef GenIntegralSet< GenSymmOper< OperatorProperties<2,true,PermutationalSymmetry::symm> >,
    IncableBFSet,
    VectorBraket<CGShell>,
    VectorBraket<CGShell>,
    QuantumNumbers<int,0> >
  DummySymmIntegral_11_11_sq;

  typedef GenIntegralSet< GenSymmOper< OperatorProperties<2,true,PermutationalSymmetry::symm> >,
    IncableBFSet,
    VectorBraket<CGF>,
    VectorBraket<CGF>,
    QuantumNumbers<int,0> >
  DummySymmIntegral_11_11_int;
};

#endif
