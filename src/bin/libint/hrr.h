
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

  /** A generic Horizontal Recurrence Relation:

  |a b) = |a+1 b-1) + AB |a b-1)

  Int is the integral class. part specifies for which particle
  the angular momentum is shifted. Function a is assumed to gain quanta,
  function b loses quanta. loc_a and loc_b specify where
  functions a and b are located (bra or ket). pos_a and pos_b
  specify which function to be used (usually pos_a and pos_b are set
  to 0 to refer to the first function for this particle in this location).

*/
  template <template <class> class I, class BFSet, int part,
  FunctionPosition loc_a, unsigned int pos_a,
  FunctionPosition loc_b, unsigned int pos_b>
  class HRR : public RecurrenceRelation {

  public:
    typedef I<BFSet> TargetType;
    typedef I<BFSet> ChildType;

    /**
      dir specifies which quantum number of a and b is shifted.
      For example, dir can be 0 (x), 1(y), or 2(z) if F is
      a Cartesian Gaussian.
      */
    HRR(const SafePtr<TargetType>&, unsigned int dir = 0);
    ~HRR();

    const unsigned int num_children() const { return num_actual_children_; };
    /// target() returns points to the i-th child
    SafePtr<TargetType> target() { return target_; };
    /// child(i) returns points i-th child
    SafePtr<ChildType> child(unsigned int i);

    const std::string cpp_function_name() {}
    const std::string cpp_source_name() {}
    const std::string cpp_header_name() {}
    std::ostream& cpp_source(std::ostream&) {}

  private:
    static const unsigned int nchild_ = 2;
    unsigned int dir_;

    SafePtr<TargetType> target_;
    SafePtr<ChildType> children_[nchild_];

    unsigned int num_actual_children_;

    void oper_checks() const;
  };

  
  
  template <template <class> class I, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    HRR<I,F,part,loc_a,pos_a,loc_b,pos_b>::HRR(const SafePtr<TargetType>& Tint, unsigned int dir) :
    target_(Tint), dir_(dir)
    {
      target_ = Tint;
      typename I<F>::AuxQuantaType aux = Tint->aux();

      typedef typename I<F>::BraType IBraType;
      typedef typename I<F>::KetType IKetType;
      IBraType* bra = new IBraType(Tint->bra());
      IKetType* ket = new IKetType(Tint->ket());

      // Zero out children pointers
      num_actual_children_ = 0;

      //
      // InBra and InKet cases have to treated explicitly since BraType and KetType don't have to match
      //
      if (loc_b == InBra) {
        // See if b-1 exists
        F sh_b(bra->member(part,pos_b));
        try {
          sh_b.dec(dir_);
        }
        catch (InvalidDecrement) {
          return;
        }
        bra->set_member(sh_b,part,pos_b);
        children_[1] = I<F>::Instance(*bra,*ket,aux);

        if (loc_a == InBra) {  // a in bra
          F sh_a(bra->member(part,pos_a));
          sh_a.inc(dir_);
          bra->set_member(sh_a,part,pos_a);
        }
        else {  // a in ket
          F sh_a(ket->member(part,pos_a));
          sh_a.inc(dir_);
          ket->set_member(sh_a,part,pos_a);
        }
        children_[0] = I<F>::Instance(*bra,*ket,aux);
        num_actual_children_ += 2;
      }
      else {
        // See if b-1 exists
        F sh_b(ket->member(part,pos_b));
        try {
          sh_b.dec(dir_);
        }
        catch (InvalidDecrement) {
          return;
        }
        ket->set_member(sh_b,part,pos_b);
        children_[1] = I<F>::Instance(*bra,*ket,aux);

        if (loc_a == InBra) {  // a in bra
          F sh_a(bra->member(part,pos_a));
          sh_a.inc(dir_);
          bra->set_member(sh_a,part,pos_a);
        }
        else {  // a in ket
          F sh_a(ket->member(part,pos_a));
          sh_a.inc(dir_);
          ket->set_member(sh_a,part,pos_a);
        }
        children_[0] = I<F>::Instance(*bra,*ket,aux);
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
    SafePtr<typename HRR<I,F,part,loc_a,pos_a,loc_b,pos_b>::ChildType>
    HRR<I,F,part,loc_a,pos_a,loc_b,pos_b>::child(unsigned int i)
    {
      assert(i>=0 && i<num_actual_children_);

      unsigned int nc=0;
      for(int c=0; c<nchild_; c++) {
        if (children_[c]) {
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
