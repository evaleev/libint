
#ifndef _libint2_src_bin_libint_dummyintegral_h_
#define _libint2_src_bin_libint_dummyintegral_h_

#include <integral.h>

namespace libint2 {
  typedef GenIntegralSet< GenMultSymm2BodyOper,
    IncableBFSet,
    DefaultTwoPBraket<CGShell>::Result,
    DefaultTwoPBraket<CGShell>::Result,
    EmptySet >
  DummySymmIntegral_11_11_sq;

  typedef GenIntegralSet< GenMultSymm2BodyOper,
    IncableBFSet,
    DefaultTwoPBraket<CGF>::Result,
    DefaultTwoPBraket<CGF>::Result,
    EmptySet >
  DummySymmIntegral_11_11_int;
};

#endif
