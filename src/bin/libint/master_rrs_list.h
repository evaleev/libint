#ifndef _libint2_src_bin_libint_masterrrslist_h_
#define _libint2_src_bin_libint_masterrrslist_h_

#include <boost/mpl/list.hpp>

#include <master_ints_list.h>
#include <hrr.h>
#include <vrr_11_twoprep_11.h>
#include <vrr_11_r12kg12_11.h>
#include <itr_11_twoprep_11.h>
#include <comp_11_r12kr12lg12_11.h>
#include <comp_11_tig12_11.h>
#include <comp_11_DivG12prime_xTx_11.h>
#include <generic_rr.h>

// master list of types and typedefs that describe capabilities of Libint2
// should not be used unless absolutely necessary.
// exception: strategy.cc and build_libint.cc, which depend on the entire typelist

namespace libint2 {

  ///////////
  // RRs
  ///////////
  typedef HRR<TwoPRep_11_11_sq,CGShell,0,InBra,0,InKet,0> HRR_ab_11_TwoPRep_11_sh;
  typedef HRR<TwoPRep_11_11_sq,CGShell,1,InBra,0,InKet,0> HRR_cd_11_TwoPRep_11_sh;
  typedef HRR<R12kG12_11_11_sq,CGShell,0,InBra,0,InKet,0> HRR_ab_11_R12kG12_11_sh;
  typedef HRR<R12kG12_11_11_sq,CGShell,1,InBra,0,InKet,0> HRR_cd_11_R12kG12_11_sh;
  typedef HRR<R12kR12lG12_11_11_sq,CGShell,0,InBra,0,InKet,0> HRR_ab_11_R12kR12lG12_11_sh;
  typedef HRR<R12kR12lG12_11_11_sq,CGShell,1,InBra,0,InKet,0> HRR_cd_11_R12kR12lG12_11_sh;
  typedef HRR<TiG12_11_11_sq,CGShell,0,InBra,0,InKet,0> HRR_ab_11_TiG12_11_sh;
  typedef HRR<TiG12_11_11_sq,CGShell,1,InBra,0,InKet,0> HRR_cd_11_TiG12_11_sh;
  typedef HRR<DivG12prime_xTx_11_11_sq,CGShell,0,InBra,0,InKet,0> HRR_ab_11_DivG12prime_xTx_sh;
  typedef HRR<DivG12prime_xTx_11_11_sq,CGShell,1,InBra,0,InKet,0> HRR_cd_11_DivG12prime_xTx_sh;

  typedef HRR<TwoPRep_11_11_int,CGF,0,InBra,0,InKet,0> HRR_ab_11_TwoPRep_11_int;
  typedef HRR<TwoPRep_11_11_int,CGF,1,InBra,0,InKet,0> HRR_cd_11_TwoPRep_11_int;
  typedef HRR<R12kG12_11_11_int,CGF,0,InBra,0,InKet,0> HRR_ab_11_R12kG12_11_int;
  typedef HRR<R12kG12_11_11_int,CGF,1,InBra,0,InKet,0> HRR_cd_11_R12kG12_11_int;
  typedef HRR<R12kR12lG12_11_11_int,CGF,0,InBra,0,InKet,0> HRR_ab_11_R12kR12lG12_11_int;
  typedef HRR<R12kR12lG12_11_11_int,CGF,1,InBra,0,InKet,0> HRR_cd_11_R12kR12lG12_11_int;
  typedef HRR<TiG12_11_11_int,CGF,0,InBra,0,InKet,0> HRR_ab_11_TiG12_11_int;
  typedef HRR<TiG12_11_11_int,CGF,1,InBra,0,InKet,0> HRR_cd_11_TiG12_11_int;
  typedef HRR<DivG12prime_xTx_11_11_int,CGF,0,InBra,0,InKet,0> HRR_ab_11_DivG12prime_xTx_int;
  typedef HRR<DivG12prime_xTx_11_11_int,CGF,1,InBra,0,InKet,0> HRR_cd_11_DivG12prime_xTx_int;

  typedef HRR<DummySymmIntegral_11_11_sq,CGShell,0,InBra,0,InKet,0> HRR_ab_11_Dummy_11_sh;
  typedef HRR<DummySymmIntegral_11_11_sq,CGShell,1,InBra,0,InKet,0> HRR_cd_11_Dummy_11_sh;
  typedef HRR<DummySymmIntegral_11_11_int,CGF,0,InBra,0,InKet,0> HRR_ab_11_Dummy_11_int;
  typedef HRR<DummySymmIntegral_11_11_int,CGF,1,InBra,0,InKet,0> HRR_cd_11_Dummy_11_int;

  typedef VRR_11_TwoPRep_11<CGShell,0,InBra> VRR_a_11_TwoPRep_11_sh;
  typedef VRR_11_TwoPRep_11<CGShell,1,InBra> VRR_c_11_TwoPRep_11_sh;
  typedef VRR_11_TwoPRep_11<CGF,0,InBra> VRR_a_11_TwoPRep_11_int;
  typedef VRR_11_TwoPRep_11<CGF,1,InBra> VRR_c_11_TwoPRep_11_int;

  typedef VRR_11_R12kG12_11<CGShell,0,InBra> VRR_a_11_R12kG12_11_sh;
  typedef VRR_11_R12kG12_11<CGShell,1,InBra> VRR_c_11_R12kG12_11_sh;
  typedef VRR_11_R12kG12_11<CGF,0,InBra> VRR_a_11_R12kG12_11_int;
  typedef VRR_11_R12kG12_11<CGF,1,InBra> VRR_c_11_R12kG12_11_int;

  typedef CR_11_R12kR12lG12_11<CGShell> CR_11_R12kR12lG12_11_sh;
  typedef CR_11_R12kR12lG12_11<CGF> CR_11_R12kR12lG12_11_int;
  
  typedef CR_11_TiG12_11<CGShell> CR_11_TiG12_11_sh;
  typedef CR_11_TiG12_11<CGF> CR_11_TiG12_11_int;

  typedef CR_11_DivG12prime_xTx_11<CGShell> CR_11_DivG12prime_xTx_11_sh;
  typedef CR_11_DivG12prime_xTx_11<CGF> CR_11_DivG12prime_xTx_11_int;

  typedef ITR_11_TwoPRep_11<GenIntegralSet_11_11,CGShell,0,InBra> ITR_a_11_TwoPRep_11_sh;
  typedef ITR_11_TwoPRep_11<GenIntegralSet_11_11,CGF,0,InBra> ITR_a_11_TwoPRep_11_int;


};

#endif // header guard
