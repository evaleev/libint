
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
  
  template <template <class> class I, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    HRR<I,F,part,loc_a,pos_a,loc_b,pos_b>::HRR(I<F>* Tint) :
    target_(Tint)
    {
      target_ = Tint;
      unsigned int m = Tint->m();

      typedef typename I<F>::BraType IBraType;
      typedef typename I<F>::KetType IKetType;
      IBraType* bra = new IBraType(Tint->bra());
      IKetType* ket = new IKetType(Tint->ket());

      // Zero out children pointers
      for(int i=0; i<nchild_; i++)
        children_[i] = 0;
      num_actual_children_ = 0;

      //
      // InBra and InKet cases have to breated explicitly since BraType and KetType don't have to match
      //
      if (loc_b == InBra) {
        // See if b-1 exists
        F sh_b(bra->member(part,pos_b));
        try {
          sh_b.dec();
        }
        catch (InvalidDecrement) {
          return;
        }
        bra->set_member(sh_b,part,pos_b);
        children_[1] = I<F>::Instance(*bra,*ket,m);

        if (loc_a == InBra) {  // a in bra
          F sh_a(bra->member(part,pos_a));
          sh_a.inc();
          bra->set_member(sh_a,part,pos_a);
        }
        else {  // a in ket
          F sh_a(ket->member(part,pos_a));
          sh_a.inc();
          ket->set_member(sh_a,part,pos_a);
        }
        children_[0] = I<F>::Instance(*bra,*ket,m);
        num_actual_children_ += 2;
      }
      else {
        // See if b-1 exists
        F sh_b(ket->member(part,pos_b));
        try {
          sh_b.dec();
        }
        catch (InvalidDecrement) {
          return;
        }
        ket->set_member(sh_b,part,pos_b);
        children_[1] = I<F>::Instance(*bra,*ket,m);

        if (loc_a == InBra) {  // a in bra
          F sh_a(bra->member(part,pos_a));
          sh_a.inc();
          bra->set_member(sh_a,part,pos_a);
        }
        else {  // a in ket
          F sh_a(ket->member(part,pos_a));
          sh_a.inc();
          ket->set_member(sh_a,part,pos_a);
        }
        children_[0] = I<F>::Instance(*bra,*ket,m);
        num_actual_children_ += 2;
      }
    }

  template <template <class> class I, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    HRR<I,F,part,loc_a,pos_a,loc_b,pos_b>::~HRR()
    {
      oper_checks();
    }

  template <template <class> class I, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    void
    HRR<I,F,part,loc_a,pos_a,loc_b,pos_b>::oper_checks() const
    {
      //
      // Here we check basic HRR applicability requirements on the integral class
      //

      // part is within the range
      typedef typename I<F>::OperatorType Oper;
      if (part < 0 || part >= Oper::Properties::np) {
        assert(false);
      }

      // can move across operator only if it's multiplicative
      if (loc_a != loc_b && !Oper::Properties::multiplicative) {
        assert(false);
      }

      // Cannot apply when a and b are the same
      if (loc_a == loc_b && pos_a == pos_b) {
        assert(false);
      }
    }
          
  template <template <class> class I, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    I<F>*
    HRR<I,F,part,loc_a,pos_a,loc_b,pos_b>::child(unsigned int i)
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

  typedef HRR<TwoPRep_11_11,CGShell,0,InBra,0,InKet,0> HRR_ab_11_TwoPRep_11_sh;
  typedef HRR<TwoPRep_11_11,CGShell,1,InBra,0,InKet,0> HRR_cd_11_TwoPRep_11_sh;
  typedef HRR<TwoPRep_11_11,CGShell,0,InKet,0,InBra,0> HRR_ba_11_TwoPRep_11_sh;
  typedef HRR<TwoPRep_11_11,CGShell,1,InKet,0,InBra,0> HRR_dc_11_TwoPRep_11_sh;


};

#endif
