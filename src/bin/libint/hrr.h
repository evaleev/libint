
#ifndef _libint2_src_bin_libint_hrr_h_
#define _libint2_src_bin_libint_hrr_h_

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <rr.h>
#include <integral.h>

using namespace std;


namespace libint2 {

  template <template <class> class I, class F, int part>
    HRR<I,F,part>::HRR(I<F>* Tint) :
    target_(Tint)
    {
      target_ = Tint;

      F sh_a(Tint->bra(0,0));
      F sh_b(Tint->ket(0,0));
      F sh_c(Tint->bra(1,0));
      F sh_d(Tint->ket(1,0));
      unsigned int m = Tint->m();

      vector<F> bra;
      vector<F> ket;

      bra.push_back(sh_a);
      bra.push_back(sh_c);
      ket.push_back(sh_b);
      ket.push_back(sh_d);

      // Zero out children pointers
      for(int i=0; i<nchild_; i++)
        children_[i] = 0;
      num_actual_children_ = 0;

      // On which particle to act
      int p_a = part;
      int p_c = (p_a == 0) ? 1 : 0;

      // See if b-1 exists
      try {
        ket[p_a].dec();
      }
      catch (InvalidDecrement) {
        return;
      }
      children_[1] = I<F>::Instance(bra[p_a],ket[p_a],bra[p_c],ket[p_c],m);
      bra[p_a].inc();
      children_[0] = I<F>::Instance(bra[p_a],ket[p_a],bra[p_c],ket[p_c],m);
      num_actual_children_ += 2;

    };

  template <template <class> class I, class F, int part>
    HRR<I,F,part>::~HRR()
    {
      if (part < 0 || part >= 2) {
        assert(false);
      }
    };

  template <template <class> class I, class F, int part>
    I<F>*
    HRR<I,F,part>::child(unsigned int i)
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

  typedef HRR<TwoERep_11_11,CGShell,0> HRR_ab_11_TwoPRep_11_sh;
  typedef HRR<TwoERep_11_11,CGShell,1> HRR_cd_11_TwoPRep_11_sh;


};

#endif
