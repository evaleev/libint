
#ifndef _libint2_src_bin_libint_itr11twoprep11_h_
#define _libint2_src_bin_libint_itr11twoprep11_h_

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <dgvertex.h>
#include <rr.h>
#include <twoprep_11_11.h>
#include <algebra.h>
#include <flop.h>
#include <prefactors.h>
#include <context.h>
#include <default_params.h>
#include <util.h>

using namespace std;


namespace libint2 {

  /** ITR (Interelectron Transfer Relation) for 2-e ERI. part specifies for which particle
  the angular momentum is raised. where specifies whether the angular momentum
  is shifted in bra or ket. Class ERI specifies which particular implementation
  of ERI to use.
  */
  template <template <typename,typename,typename> class ERI, class BFSet, int part, FunctionPosition where>
    class ITR_11_TwoPRep_11 : public RecurrenceRelation
    {

  public:
    typedef RecurrenceRelation ParentType;
    typedef BFSet BasisFunctionType;
    typedef ITR_11_TwoPRep_11 ThisType;
    typedef ERI<BFSet,TwoPRep,mType> TargetType;
    typedef TargetType ChildType;
    /// The type of expressions in which RecurrenceRelations result.
    typedef RecurrenceRelation::ExprType ExprType;

    /** Use Instance() to obtain an instance of RR. This function is provided to avoid
        issues with getting a SafePtr from constructor (as needed for registry to work).

        dir specifies which quantum number of a and b is shifted.
        For example, dir can be 0 (x), 1(y), or 2(z) if F is
        a Cartesian Gaussian.
    */
    static SafePtr<ThisType> Instance(const SafePtr<TargetType>&, unsigned int dir = 0);
    ~ITR_11_TwoPRep_11() { assert(part == 0 || part == 1); }

    /// Implementation of RecurrenceRelation::num_children()
    const unsigned int num_children() const { return children_.size(); };
    /// Implementation of RecurrenceRelation::rr_target()
    SafePtr<DGVertex> rr_target() const { return static_pointer_cast<DGVertex,TargetType>(target_); }
    /// Implementation of RecurrenceRelation::rr_child()
    SafePtr<DGVertex> rr_child(unsigned int i) const { return static_pointer_cast<DGVertex,ChildType>(children_.at(i)); }
    /// Implementation of RecurrenceRelation::is_simple()
    bool is_simple() const {
      return TrivialBFSet<BFSet>::result;
    }

  private:
    /**
      dir specifies which quantum number is incremented.
      For example, dir can be 0 (x), 1(y), or 2(z) if BFSet is
      a Cartesian Gaussian.
     */
    ITR_11_TwoPRep_11(const SafePtr<TargetType>&, unsigned int dir);

    static const unsigned int max_nchildren_ = 4;
    unsigned int dir_;
    SafePtr<TargetType> target_;
    std::vector< SafePtr<ChildType> > children_;
    const SafePtr<ChildType>& make_child(const BFSet& A, const BFSet& B, const BFSet& C, const BFSet& D, unsigned int m) {
      const SafePtr<ChildType>& i = ChildType::Instance(A,B,C,D,m);
      children_.push_back(i);
      return *(children_.end()-1);
    }

    std::string generate_label() const
    {
      typedef typename TargetType::AuxIndexType mType;
      static SafePtr<mType> aux0(new mType(0u));
      ostringstream os;
      // ITR recurrence relation code is independent of m (it never appears anywhere in equations), hence
      // to avoid generating identical code make sure that the (unique) label does not contain m
      os << "ITR Part" << part << " " << to_string(where) << genintegralset_label(target_->bra(),target_->ket(),aux0,target_->oper());
      return os.str();
    }

  };

  template <template <typename,typename,typename> class ERI, class F, int part, FunctionPosition where>
    SafePtr< ITR_11_TwoPRep_11<ERI,F,part,where> >
    ITR_11_TwoPRep_11<ERI,F,part,where>::Instance(const SafePtr<TargetType>& Tint,
                                                  unsigned int dir)
    {
      SafePtr<ThisType> this_ptr(new ThisType(Tint,dir));
      // Do post-construction duties
      if (this_ptr->num_children() != 0) {
        this_ptr->register_with_rrstack<ThisType>();
        return this_ptr;
      }
      return SafePtr<ThisType>();
    }

  template <template <typename,typename,typename> class ERI, class F, int part, FunctionPosition where>
    ITR_11_TwoPRep_11<ERI,F,part,where>::ITR_11_TwoPRep_11(const SafePtr<TargetType>& Tint,
                                                           unsigned int dir) :
    target_(Tint), dir_(dir)
    {
      /// TODO debug this RR, it seems broken (MPQC libint2 tests fail in HF)
      throw ProgrammingError("ITR_11_TwoPRep_11 is broken");

      /// InKet
      if (where == InKet)
        throw ProgrammingError("ITR_11_TwoPRep_11<ERI,F,part,where>::ITR_11_TwoPRep_11() -- where=InKet not implemented yet");

      {
        F a(Tint->bra(0,0));
        F b(Tint->ket(0,0));
        F c(Tint->bra(1,0));
        F d(Tint->ket(1,0));
        if (a.contracted() ||
            b.contracted() ||
            c.contracted() ||
            d.contracted())
          return;
      }

      children_.reserve(max_nchildren_);
      using namespace libint2::algebra;
      using namespace libint2::prefactor;
      const unsigned int m = Tint->aux()->elem(0);
      const F& _1 = unit<F>(dir);

      // B and D must be s functions -- general RR not implemented yet
      {
        F b(Tint->ket(0,0) - _1);
        F d(Tint->ket(1,0) - _1);
        if (exists(b) || exists(d))
          return;
      }

      // Build on A
      if (part == 0 && where == InBra) {
        F a(Tint->bra(0,0) - _1);
        if (!exists(a)) return;
        F b(Tint->ket(0,0));
        F c(Tint->bra(1,0));
        F d(Tint->ket(1,0));

        const SafePtr<ChildType>& ABCD_m = make_child(a,b,c,d,m);
        if (is_simple()) { expr_ = Vector("TwoPRepITR_pfac0_0")[dir] * ABCD_m;  nflops_+=1; }

        const F& cp1 = c + _1;
        const SafePtr<ChildType>& ABCp1D_m = make_child(a,b,cp1,d,m);
        if (is_simple()) { expr_ += Scalar("TwoPRepITR_pfac1_0") * ABCp1D_m;  nflops_+=2; }

        const F& am1 = a - _1;
        if (exists(am1)) {
          const SafePtr<ChildType>& Am1BCD_m = make_child(am1,b,c,d,m);
          if (is_simple()) { expr_ += Vector(a)[dir] * Scalar("oo2z") * Am1BCD_m;  nflops_+=3; }
        }
        const F& cm1 = c - _1;
        if (exists(cm1)) {
          const SafePtr<ChildType>& ABCm1D_m = make_child(a,b,cm1,d,m);
          if (is_simple()) { expr_ += Vector(c)[dir] * Scalar("oo2z") * ABCm1D_m;  nflops_+=3; }
        }
        return;
      }
      // Build on C
      if (part == 1 && where == InBra) {
        F a(Tint->bra(0,0));
        F b(Tint->ket(0,0));
        F c(Tint->bra(1,0) - _1);
        if (!exists(c)) return;
        F d(Tint->ket(1,0));

        const SafePtr<ChildType>& ABCD_m = make_child(a,b,c,d,m);
        if (is_simple()) { expr_ = Vector("TwoPRepITR_pfac0_1")[dir] * ABCD_m;  nflops_+=1; }

        const F& ap1 = a + _1;
        const SafePtr<ChildType>& Ap1BCD_m = make_child(ap1,b,c,d,m);
        if (is_simple()) { expr_ += Scalar("TwoPRepITR_pfac1_1") * Ap1BCD_m;  nflops_+=2; }

        const F& cm1 = c - _1;
        if (exists(cm1)) {
          const SafePtr<ChildType>& ABCm1D_m = make_child(a,b,cm1,d,m);
          if (is_simple()) { expr_ += Vector(c)[dir] * Scalar("oo2e") * ABCm1D_m;  nflops_+=3; }
        }
        const F& am1 = a - _1;
        if (exists(am1)) {
          const SafePtr<ChildType>& Am1BCD_m = make_child(am1,b,c,d,m);
          if (is_simple()) { expr_ += Vector(a)[dir] * Scalar("oo2e") * Am1BCD_m;  nflops_+=3; }
        }
        return;
      }
    }

};

#endif
