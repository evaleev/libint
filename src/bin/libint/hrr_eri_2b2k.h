
#ifndef _libint2_src_bin_libint_hrreri2b2k_h_
#define _libint2_src_bin_libint_hrreri2b2k_h_

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <rr.h>

using namespace std;


namespace libint2 {

  template <template <class> class I, class F, int part>
    HRR_ERI_2b2k<I,F,part>::HRR_ERI_2b2k(I<F>* Tint) :
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

      // On which particle to act
      int p_a = part;
      int p_c = (p_a == 0) ? 1 : 0;

      // See if b-1 exists
      try {
        ket[p_a][0].dec();
      }
      catch (InvalidDecrement) {
        return;
      }
      children_[1] = I<F>::Instance(bra,ket,m);
      bra[p_a][0].inc();
      children_[0] = I<F>::Instance(bra,ket,m);
      num_actual_children_ += 2;

    };

  template <template <class> class I, class F, int part>
    HRR_ERI_2b2k<I,F,part>::~HRR_ERI_2b2k()
    {
      if (part < 0 || part >= 2) {
        assert(false);
      }
    };

  template <template <class> class I, class F, int part>
    I<F>*
    HRR_ERI_2b2k<I,F,part>::child(unsigned int i)
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

  typedef HRR_ERI_2b2k<TwoERep_2b2k,CGShell,0> HRR_ab_ERI_2b2k_shell;
  typedef HRR_ERI_2b2k<TwoERep_2b2k,CGShell,1> HRR_cd_ERI_2b2k_shell;


};

#endif
