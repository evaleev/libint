
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

  template <template <class> class ERI, class F, int part, bool use_bra>
    VRR_ERI_2b2k<ERI,F,part,use_bra>::VRR_ERI_2b2k(ERI<F>* Tint) :
    target_(Tint)
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
      for(int i=0; i<nchild_; i++)
        children_[i] = 0;
      num_actual_children_ = 0;

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
      children_[0] = ERI<F>::Instance(bra,ket,m);
      children_[1] = ERI<F>::Instance(bra,ket,m+1);
      num_actual_children_ += 2;
      // See if a-2 exists
      bool a_minus_2_exists = true;
      try {
        braket[p_a][0].dec();
      }
      catch (InvalidDecrement) {
        a_minus_2_exists = false;
      }
      if (a_minus_2_exists) {
        children_[2] = ERI<F>::Instance(bra,ket,m);
        children_[3] = ERI<F>::Instance(bra,ket,m+1);
        braket[p_a][0].inc();
        num_actual_children_ += 2;
      }

      try {
        braket[p_c][0].dec();
      }
      catch (InvalidDecrement) {
        return;
      }
      children_[4] = ERI<F>::Instance(bra,ket,m+1);
      num_actual_children_ += 1;

    };

  template <template <class> class ERI, class F, int part, bool use_bra>
    VRR_ERI_2b2k<ERI,F,part,use_bra>::~VRR_ERI_2b2k()
    {
      if (part < 0 || part >= 2) {
        assert(false);
      }
    };

  template <template <class> class ERI, class F, int part, bool use_bra>
    ERI<F>*
    VRR_ERI_2b2k<ERI,F,part,use_bra>::child(unsigned int i)
    {
      assert(i>=0 && i<num_actual_children_);

      unsigned int nc=0;
      for(int c=0; c<nchild_; c++) {
        if (children_[c] != 0) {
          if (nc == i)
            return children_[c];
          nc++;
        }
      }
    };

  typedef VRR_ERI_2b2k<TwoERep_2b2k,CGShell,0,true> VRR_a_ERI_2b2k_shell;
  typedef VRR_ERI_2b2k<TwoERep_2b2k,CGShell,1,true> VRR_c_ERI_2b2k_shell;
  typedef VRR_ERI_2b2k<TwoERep_2b2k,CGShell,0,false> VRR_b_ERI_2b2k_shell;
  typedef VRR_ERI_2b2k<TwoERep_2b2k,CGShell,1,false> VRR_d_ERI_2b2k_shell;


};

#endif
