
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

  template <template <class> class I, class F, int part, bool to_bra>
    HRR<I,F,part,to_bra>::HRR(I<F>* Tint) :
    target_(Tint)
    {
      target_ = Tint;
      unsigned int m = Tint->m();

      typedef typename I<F>::BraType IBraType;
      typedef typename I<F>::KetType IKetType;
      IBraType* bra = Tint->bra();
      IKetType* ket = Tint->ket();

      // Zero out children pointers
      for(int i=0; i<nchild_; i++)
        children_[i] = 0;
      num_actual_children_ = 0;

      // On which particle to act
      int p_a = part;

      if (to_bra) {
        // See if b-1 exists
        F sh_b(ket->member(p_a,0));
        try {
          sh_b.dec();
        }
        catch (InvalidDecrement) {
          return;
        }
        ket->set_member(sh_b,p_a,0);
        children_[1] = I<F>::Instance(*bra,*ket,m);

        F sh_a(bra->member(p_a,0));
        sh_a.inc();
        bra->set_member(sh_a,p_a,0);
        children_[0] = I<F>::Instance(*bra,*ket,m);
        num_actual_children_ += 2;
      }
      // to_ket
      else {
        // See if a-1 exists
        F sh_a(bra->member(p_a,0));
        try {
          sh_a.dec();
        }
        catch (InvalidDecrement) {
          return;
        }
        bra->set_member(sh_a,p_a,0);
        children_[1] = I<F>::Instance(*bra,*ket,m);

        F sh_b(ket->member(p_a,0));
        sh_b.inc();
        ket->set_member(sh_b,p_a,0);
        children_[0] = I<F>::Instance(*bra,*ket,m);
        num_actual_children_ += 2;
      }

    };

  template <template <class> class I, class F, int part, bool to_bra>
    HRR<I,F,part,to_bra>::~HRR()
    {
      if (part < 0 || part >= 2) {
        assert(false);
      }
    };

  template <template <class> class I, class F, int part, bool to_bra>
    I<F>*
    HRR<I,F,part,to_bra>::child(unsigned int i)
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

  typedef HRR<TwoPRep_11_11,CGShell,0,true> HRR_ab_11_TwoPRep_11_sh;
  typedef HRR<TwoPRep_11_11,CGShell,1,true> HRR_cd_11_TwoPRep_11_sh;
  typedef HRR<TwoPRep_11_11,CGShell,0,false> HRR_ba_11_TwoPRep_11_sh;
  typedef HRR<TwoPRep_11_11,CGShell,1,false> HRR_dc_11_TwoPRep_11_sh;


};

#endif
