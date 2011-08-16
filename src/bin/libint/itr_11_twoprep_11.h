
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
    virtual ~ITR_11_TwoPRep_11() { assert(part == 0 || part == 1); }

#if 1
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
#endif

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

#if LIBINT_ENABLE_GENERIC_CODE
    /// Implementation of RecurrenceRelation::has_generic()
    bool has_generic(const SafePtr<CompilationParameters>& cparams) const;
    /// Implementation of RecurrenceRelation::generic_header()
    std::string generic_header() const;
    /// Implementation of RecurrenceRelation::generic_instance()
    std::string generic_instance(const SafePtr<CodeContext>& context, const SafePtr<CodeSymbols>& args) const;
#endif
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
      // only works for primitive integrals
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

      // not yet implemented for derivative integrals
      {
        const OriginDerivative dA = Tint->bra(0,0).deriv();
        const OriginDerivative dB = Tint->ket(0,0).deriv();
        const OriginDerivative dC = Tint->bra(1,0).deriv();
        const OriginDerivative dD = Tint->ket(1,0).deriv();
        const bool deriv = dA.zero() == false ||
            dB.zero() == false ||
            dC.zero() == false ||
            dD.zero() == false;
        if (deriv)
          return;
      }

      children_.reserve(max_nchildren_);
      using namespace libint2::algebra;
      using namespace libint2::prefactor;
      const unsigned int m = Tint->aux()->elem(0);
      const F& _1 = unit<F>(dir);

      // Build on A
      if (part == 0 && where == InBra) {
        F a(Tint->bra(0,0) - _1);
        if (!exists(a)) return;
        F b(Tint->ket(0,0));
        F c(Tint->bra(1,0));
        F d(Tint->ket(1,0));

        const SafePtr<ChildType>& ABCD_m = make_child(a,b,c,d,m);
        if (is_simple()) { expr_ = Vector("TwoPRepITR_pfac0_0_0")[dir] * ABCD_m;  nflops_+=1; }

        const F& cp1 = c + _1;
        const SafePtr<ChildType>& ABCp1D_m = make_child(a,b,cp1,d,m);
        if (is_simple()) { expr_ -= Scalar("eoz") * ABCp1D_m;  nflops_+=2; }

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
        if (is_simple()) { expr_ = Vector("TwoPRepITR_pfac0_1_0")[dir] * ABCD_m;  nflops_+=1; }

        const F& ap1 = a + _1;
        const SafePtr<ChildType>& Ap1BCD_m = make_child(ap1,b,c,d,m);
        if (is_simple()) { expr_ -= Scalar("zoe") * Ap1BCD_m;  nflops_+=2; }

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
      // Build on B
      if (part == 0 && where == InKet) {
        F a(Tint->bra(0,0));
        F b(Tint->ket(0,0) - _1);
        if (!exists(b)) return;
        F c(Tint->bra(1,0));
        F d(Tint->ket(1,0));

        const SafePtr<ChildType>& ABCD_m = make_child(a,b,c,d,m);
        if (is_simple()) { expr_ = Vector("TwoPRepITR_pfac0_0_1")[dir] * ABCD_m;  nflops_+=1; }

        const F& dp1 = d + _1;
        const SafePtr<ChildType>& ABCDp1_m = make_child(a,b,c,dp1,m);
        if (is_simple()) { expr_ -= Scalar("eoz") * ABCDp1_m;  nflops_+=2; }

        const F& bm1 = b - _1;
        if (exists(bm1)) {
          const SafePtr<ChildType>& ABm1CD_m = make_child(a,bm1,c,d,m);
          if (is_simple()) { expr_ += Vector(b)[dir] * Scalar("oo2z") * ABm1CD_m;  nflops_+=3; }
        }
        const F& dm1 = d - _1;
        if (exists(dm1)) {
          const SafePtr<ChildType>& ABCDm1_m = make_child(a,b,c,dm1,m);
          if (is_simple()) { expr_ += Vector(d)[dir] * Scalar("oo2z") * ABCDm1_m;  nflops_+=3; }
        }
        return;
      }
      // Build on D
      if (part == 1 && where == InBra) {
        F a(Tint->bra(0,0));
        F b(Tint->ket(0,0));
        F c(Tint->bra(1,0));
        F d(Tint->ket(1,0) - _1);
        if (!exists(d)) return;

        const SafePtr<ChildType>& ABCD_m = make_child(a,b,c,d,m);
        if (is_simple()) { expr_ = Vector("TwoPRepITR_pfac0_1_1")[dir] * ABCD_m;  nflops_+=1; }

        const F& bp1 = b + _1;
        const SafePtr<ChildType>& ABp1CD_m = make_child(a,bp1,c,d,m);
        if (is_simple()) { expr_ -= Scalar("zoe") * ABp1CD_m;  nflops_+=2; }

        const F& dm1 = d - _1;
        if (exists(dm1)) {
          const SafePtr<ChildType>& ABCDm1_m = make_child(a,b,c,dm1,m);
          if (is_simple()) { expr_ += Vector(d)[dir] * Scalar("oo2e") * ABCDm1_m;  nflops_+=3; }
        }
        const F& bm1 = b - _1;
        if (exists(bm1)) {
          const SafePtr<ChildType>& ABm1CD_m = make_child(a,bm1,c,d,m);
          if (is_simple()) { expr_ += Vector(b)[dir] * Scalar("oo2e") * ABm1CD_m;  nflops_+=3; }
        }
        return;
      }

    }

#if LIBINT_ENABLE_GENERIC_CODE
  template <template <typename,typename,typename> class ERI, class F, int part, FunctionPosition where>
  bool
    ITR_11_TwoPRep_11<ERI,F,part,where>::has_generic(const SafePtr<CompilationParameters>& cparams) const
    {
      if (TrivialBFSet<F>::result)
        return false;

      F sh_a(target_->bra(0,0));
      F sh_b(target_->ket(0,0));
      F sh_c(target_->bra(1,0));
      F sh_d(target_->ket(1,0));
      const unsigned int max_opt_am = cparams->max_am_opt();
      // generic code works for a0c0 of 0a0c classes where am(a) > 1 and am(c) > 1
      // to generate optimized code for xxxx integral need to generate specialized code for up to (x+x)0(x+x)0 integrals
      if (sh_b.zero() && sh_d.zero() &&
          (sh_a.norm() > std::max(2*max_opt_am,1u) ||
           sh_c.norm() > std::max(2*max_opt_am,1u)
          ) &&
          (sh_a.norm() > 1u && sh_c.norm() > 1u)
         )
        return true;
      if (sh_a.zero() && sh_c.zero() &&
          (sh_b.norm() > std::max(2*max_opt_am,1u) ||
           sh_d.norm() > std::max(2*max_opt_am,1u)
          ) &&
          (sh_b.norm() > 1u && sh_d.norm() > 1u)
         )
        return true;
      return false;
    }

  template <template <typename,typename,typename> class ERI, class F, int part, FunctionPosition where>
  std::string
    ITR_11_TwoPRep_11<ERI,F,part,where>::generic_header() const
    {
      F sh_a(target_->bra(0,0));
      F sh_b(target_->ket(0,0));
      F sh_c(target_->bra(1,0));
      F sh_d(target_->ket(1,0));
      const bool xsxs = sh_b.zero() && sh_d.zero();
      const bool sxsx = sh_a.zero() && sh_c.zero();

      const OriginDerivative dA = target_->bra(0,0).deriv();
      const OriginDerivative dB = target_->ket(0,0).deriv();
      const OriginDerivative dC = target_->bra(1,0).deriv();
      const OriginDerivative dD = target_->ket(1,0).deriv();
      const bool deriv = dA.zero() == false ||
          dB.zero() == false ||
          dC.zero() == false ||
          dD.zero() == false;
      assert(deriv == false);

      if (deriv == false) {
        return std::string("ITR_xs_xs.h");
      }
      else {
        if (xsxs) return std::string("ITR_xs_xs_deriv.h");
        if (sxsx) return std::string("ITR_sx_sx_deriv.h");
      }
      abort(); // unreachable
    }

  template <template <typename,typename,typename> class ERI, class F, int part, FunctionPosition where>
  std::string
    ITR_11_TwoPRep_11<ERI,F,part,where>::generic_instance(const SafePtr<CodeContext>& context, const SafePtr<CodeSymbols>& args) const
    {
      std::ostringstream oss;
      F sh_a(target_->bra(0,0));
      F sh_b(target_->ket(0,0));
      F sh_c(target_->bra(1,0));
      F sh_d(target_->ket(1,0));
      const bool xsxs = sh_b.zero() && sh_d.zero();
      const bool sxsx = sh_a.zero() && sh_c.zero();

      const OriginDerivative dA = target_->bra(0,0).deriv();
      const OriginDerivative dB = target_->ket(0,0).deriv();
      const OriginDerivative dC = target_->bra(1,0).deriv();
      const OriginDerivative dD = target_->ket(1,0).deriv();
      const bool deriv = dA.zero() == false ||
          dB.zero() == false ||
          dC.zero() == false ||
          dD.zero() == false;
      assert(deriv == false);

      oss << "using namespace libint2;" << endl;

      if (deriv == false) { // for regular integrals I know exactly how many prerequisites I need
        if(xsxs) {
          oss << "libint2::ITR_xs_xs<" << part << "," << sh_a.norm() << "," << sh_c.norm() << ",true,";
        }
        if (sxsx) {
          oss << "libint2::ITR_xs_xs<" << part << "," << sh_b.norm() << "," << sh_d.norm() << ",false,";
        }
        oss << ((context->cparams()->max_vector_length() == 1) ? "false" : "true");
        oss << ">::compute(inteval";

        const unsigned int nargs = args->n();
        for(unsigned int a=0; a<nargs; a++) {
          oss << "," << args->symbol(a);
        }
        oss << ");";
      }

      return oss.str();
    }
#endif // #if !LIBINT_ENABLE_GENERIC_CODE

};

#endif
