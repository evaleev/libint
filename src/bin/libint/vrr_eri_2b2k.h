
#ifndef _libint2_src_bin_libint_vrreri2b2k_h_
#define _libint2_src_bin_libint_vrreri2b2k_h_

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <rr.h>

using namespace std;


namespace libint2 {

  template <class F, int part, bool use_bra>
    VRR_ERI_2b2k<F,part,use_bra>::VRR_ERI_2b2k(const TwoERep_2b2k<F>* Tint)
    {
      target_ = Tint;

      F sh_a = Tint->bra(0,0);
      F sh_b = Tint->ket(0,0);
      F sh_c = Tint->bra(1,0);
      F sh_d = Tint->ket(1,0);
      unsigned int m = Tint->m();

      vector<F> bra[2];
      vector<F> ket[2];

      bra[0].push_back(sh_a);
      bra[1].push_back(sh_c);
      ket[0].push_back(sh_b);
      ket[1].push_back(sh_d);

      // Zero out children pointers
      for(int i=0; i<5; i++)
        children_[i] = 0;


      // Use indirection to choose bra or ket
      vector<F>* braket;
      if (use_bra) {
        braket = bra;
      }
      else {
        braket = ket;
      }
      // On which particle to act
      int p_a = part;
      int p_c = (p_a == 0) ? 1 : 0;

      // See if a-1 exists
      try {
        braket[p_a][0].dec();
      }
      catch (InvalidDecrement) {
        return;
      }
      children_[0] = new TwoERep_2b2k<F>(bra,ket,m);
      children_[1] = new TwoERep_2b2k<F>(bra,ket,m+1);
      // See if a-2 exists
      bool a_minus_2_exists = true;
      try {
        braket[p_a][0].dec();
      }
      catch (InvalidDecrement) {
        a_minus_2_exists = false;
      }
      if (a_minus_2_exists) {
        children_[2] = new TwoERep_2b2k<F>(bra,ket,m);
        children_[3] = new TwoERep_2b2k<F>(bra,ket,m+1);
        braket[p_a][0].inc();
      }

      try {
        braket[p_c][0].dec();
      }
      catch (InvalidDecrement) {
        return;
      }
      children_[4] = new TwoERep_2b2k<F>(bra,ket,m);

    };

  template <class F, int part, bool use_bra>
    VRR_ERI_2b2k<F,part,use_bra>::~VRR_ERI_2b2k()
    {
      if (part < 0 || part >= 2) {
        assert(false);
      }
    };

  typedef VRR_ERI_2b2k<CGShell,0,true> VRR_a_ERI_2b2k_shell;
  typedef VRR_ERI_2b2k<CGShell,1,true> VRR_c_ERI_2b2k_shell;
  typedef VRR_ERI_2b2k<CGShell,0,false> VRR_b_ERI_2b2k_shell;
  typedef VRR_ERI_2b2k<CGShell,1,false> VRR_d_ERI_2b2k_shell;


};

#endif
