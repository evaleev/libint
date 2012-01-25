#ifndef _libint2_src_bin_libint_masterintslist_h_
#define _libint2_src_bin_libint_masterintslist_h_

#include <boost/mpl/list.hpp>

#include <bfset.h>
#include <oper.h>
#if LIBINT_SUPPORT_ONEBODYINTS
#  include <integral_1_1.h>
#endif
#include <integral_11_11.h>

namespace libint2 {

  //////////////////////////
  // one-electron integrals
  //////////////////////////
#if LIBINT_SUPPORT_ONEBODYINTS
  typedef GenIntegralSet_1_1<CGShell,OnePSep,EmptySet> OnePSep_1_1_sq;
  typedef GenIntegralSet_1_1<CGF,OnePSep,EmptySet> OnePSep_1_1_int;
  typedef GenIntegralSet_1_1<CGShell,OnePNonSep,mType> OnePNonSep_1_1_sq;
  typedef GenIntegralSet_1_1<CGF,OnePNonSep,mType> OnePNonSep_1_1_int;
#endif

  //////////////////////////
  // two-electron integrals
  //////////////////////////
  typedef GenIntegralSet_11_11<CGShell,TwoPRep,mType> TwoPRep_11_11_sq;
  typedef GenIntegralSet_11_11<CGF,TwoPRep,mType> TwoPRep_11_11_int;
  typedef GenIntegralSet_11_11<CGShell,R12kG12,mType> R12kG12_11_11_sq;
  typedef GenIntegralSet_11_11<CGF,R12kG12,mType> R12kG12_11_11_int;
  typedef GenIntegralSet_11_11<CGShell,R12kR12lG12,EmptySet> R12kR12lG12_11_11_sq;
  typedef GenIntegralSet_11_11<CGF,R12kR12lG12,EmptySet> R12kR12lG12_11_11_int;
  typedef GenIntegralSet_11_11<CGShell,TiG12,mType> TiG12_11_11_sq;
  typedef GenIntegralSet_11_11<CGF,TiG12,mType> TiG12_11_11_int;
  typedef GenIntegralSet_11_11<CGShell,G12TiG12,mType> G12TiG12_11_11_sq;
  typedef GenIntegralSet_11_11<CGF,G12TiG12,mType> G12TiG12_11_11_int;
  typedef GenIntegralSet_11_11<CGShell,DivG12prime_xTx,mType> DivG12prime_xTx_11_11_sq;
  typedef GenIntegralSet_11_11<CGF,DivG12prime_xTx,mType> DivG12prime_xTx_11_11_int;
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

  /** All known types go into this typelist
      Every type must have a corresponding instantiation of MasterStrategy in strategy.cc
    */
  typedef mpl::list<
#if LIBINT_SUPPORT_ONEBODYINTS
  OnePSep_1_1_sq,
  OnePSep_1_1_int,
  OnePNonSep_1_1_sq,
  OnePNonSep_1_1_int,
#endif
  TwoPRep_11_11_sq,
  TwoPRep_11_11_int,
  R12kG12_11_11_sq,
  R12kG12_11_11_int,
  R12kR12lG12_11_11_sq,
  R12kR12lG12_11_11_int,
  TiG12_11_11_sq,
  TiG12_11_11_int,
  G12TiG12_11_11_sq,
  G12TiG12_11_11_int,
  DivG12prime_xTx_11_11_sq,
  DivG12prime_xTx_11_11_int,
  DummySymmIntegral_11_11_sq,
  DummySymmIntegral_11_11_int
  > MasterIntegralTypeList;

};

#endif // header guard
