
#ifndef _libint2_src_bin_libint_vrr11twoprep11_h_
#define _libint2_src_bin_libint_vrr11twoprep11_h_

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <rr.h>
#include <integral.h>

using namespace std;


namespace libint2 {

  /** VRR Recurrence Relation for 2-e ERI. part specifies for which particle
  the angular momentum is raised. bool bra specifies whether the angular momentum
  is raised in bra (true) or ket (false). Class ERI specifies which particular implementation
  of ERI to use.
  */
  template <template <class> class ERI, class BFSet, int part, FunctionPosition where>
  class VRR_11_TwoPRep_11 : public RecurrenceRelation {
  public:
    typedef ERI<BFSet> TargetType;
    typedef ERI<BFSet> ChildType;

    /**
      dir specifies which quantum number is incremented.
      For example, dir can be 0 (x), 1(y), or 2(z) if BFSet is
      a Cartesian Gaussian.
     */
    VRR_11_TwoPRep_11(const SafePtr<TargetType>&, unsigned int dir = 0);
    ~VRR_11_TwoPRep_11();

    const unsigned int num_children() const { return num_actual_children_; };
    /// target() returns pointer to the i-th child
    SafePtr<TargetType> target() const { return target_; };
    /// child(i) returns pointer to the i-th child
    SafePtr<ChildType> child(unsigned int i) const;
    /// Implementation of RecurrenceRelation's target()
    SafePtr<DGVertex> rr_target() const { return static_pointer_cast<DGVertex,TargetType>(target()); }
    /// Implementation of RecurrenceRelation's child()
    SafePtr<DGVertex> rr_child(unsigned int i) const { return static_pointer_cast<DGVertex,ChildType>(child(i)); }
    
    const std::string cpp_function_name() {}
    const std::string cpp_source_name() {}
    const std::string cpp_header_name() {}
    std::ostream& cpp_source(std::ostream&) {}

  private:
    static const unsigned int nchild_ = 5;
    unsigned int dir_;

    SafePtr<TargetType> target_;
    SafePtr<ChildType> children_[nchild_];

    unsigned int num_actual_children_;

  };

  
  template <template <class> class ERI, class F, int part, FunctionPosition where>
    VRR_11_TwoPRep_11<ERI,F,part,where>::VRR_11_TwoPRep_11(const SafePtr<ERI<F> >& Tint,
                                                           unsigned int dir) :
    target_(Tint), dir_(dir)
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
      num_actual_children_ = 0;

      // Use indirection to choose bra or ket
      vector<F>* bra_ref = &bra;
      vector<F>* ket_ref = &ket;
      if (where == InKet) {
        bra_ref = &ket;
        ket_ref = &bra;
      }
      // On which particle to act
      int p_a = part;
      int p_c = (p_a == 0) ? 1 : 0;

      // See if a-1 exists
      try {
        bra_ref->operator[](p_a).dec(dir);
      }
      catch (InvalidDecrement) {
        return;
      }
      children_[0] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m);
      children_[1] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m+1);
      num_actual_children_ += 2;
      // See if a-2 exists
      bool a_minus_2_exists = true;
      try {
        bra_ref->operator[](p_a).dec(dir);
      }
      catch (InvalidDecrement) {
        a_minus_2_exists = false;
      }
      if (a_minus_2_exists) {
        children_[2] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m);
        children_[3] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m+1);
        bra_ref->operator[](p_a).inc(dir);
        num_actual_children_ += 2;
      }

      try {
        bra_ref->operator[](p_c).dec(dir);
      }
      catch (InvalidDecrement) {
        return;
      }
      children_[4] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m+1);
      num_actual_children_ += 1;

    };

  template <template <class> class ERI, class F, int part, FunctionPosition where>
    VRR_11_TwoPRep_11<ERI,F,part,where>::~VRR_11_TwoPRep_11()
    {
      if (part < 0 || part >= 2) {
        assert(false);
      }
    };

  template <template <class> class ERI, class F, int part, FunctionPosition where>
    SafePtr< ERI<F> >
    VRR_11_TwoPRep_11<ERI,F,part,where>::child(unsigned int i) const
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

  typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGShell,0,InBra> VRR_a_11_TwoPRep_11_sh;
  typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGShell,1,InBra> VRR_c_11_TwoPRep_11_sh;
  typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGShell,0,InKet> VRR_b_11_TwoPRep_11_sh;
  typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGShell,1,InKet> VRR_d_11_TwoPRep_11_sh;


};

#endif
