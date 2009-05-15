
#ifndef _libint2_src_bin_libint_vrr11twoprep11_h_
#define _libint2_src_bin_libint_vrr11twoprep11_h_

#include <generic_rr.h>
#include <twoprep_11_11.h>

using namespace std;

namespace libint2 {

  /** VRR Recurrence Relation for 2-e ERI. part specifies for which particle
  the angular momentum is raised. where specifies whether angular momentus is decreased in bra or ket.
  */
    template <class BFSet, int part, FunctionPosition where>
      class VRR_11_TwoPRep_11 : public GenericRecurrenceRelation< VRR_11_TwoPRep_11<BFSet,part,where>,
                                                                  BFSet,
                                                                  GenIntegralSet_11_11<BFSet,TwoPRep,mType> >
    {
    public:
      typedef VRR_11_TwoPRep_11 ThisType;
      typedef BFSet BasisFunctionType;
      typedef GenIntegralSet_11_11<BFSet,TwoPRep,mType> TargetType;
      typedef GenericRecurrenceRelation<ThisType,BFSet,TargetType> ParentType;
      friend class GenericRecurrenceRelation<ThisType,BFSet,TargetType>;
      static const unsigned int max_nchildren = 8;

      using ParentType::Instance;
    private:
      using RecurrenceRelation::expr_;
      using RecurrenceRelation::nflops_;
      using ParentType::target_;
      using ParentType::is_simple;

      /// Constructor is private, used by ParentType::Instance that mainains registry of these objects
      VRR_11_TwoPRep_11(const SafePtr<TargetType>&, unsigned int dir);
      /// Default directionality
      static bool directional() { return ParentType::default_directional(); }

      static std::string descr() { return "OSVRR"; }
      /** Re-Implementation of GenericRecurrenceRelation::generate_label():
          TwoPRep VRR recurrence relations codes are independent of m (it never appears anywhere in equations), hence
          to avoid generating identical code make sure that the (unique) label does not contain m. */
      std::string generate_label() const
      {
        typedef typename TargetType::AuxIndexType mType;
        static SafePtr<mType> aux0(new mType(0u));
        ostringstream os;
        os << descr() << "P" << part << to_string(where)
           << genintegralset_label(target_->bra(),target_->ket(),aux0,target_->oper());
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

  template <class F, int part, FunctionPosition where>
    VRR_11_TwoPRep_11<F,part,where>::VRR_11_TwoPRep_11(const SafePtr< TargetType >& Tint,
                                                       unsigned int dir) :
    ParentType(Tint,dir)
    {
      using namespace libint2::algebra;
      using namespace libint2::prefactor;
      using namespace libint2::braket;
      const unsigned int m = Tint->aux()->elem(0);
      const F& _1 = unit<F>(dir);

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

      typedef TargetType ChildType;
      ChildFactory<ThisType,ChildType> factory(this);

      // Build on A or B
      if (part == 0 && where == InBra) {
        F a(Tint->bra(0,0) - _1);
        if (!exists(a)) return;
        F b(Tint->ket(0,0));
        F c(Tint->bra(1,0));
        F d(Tint->ket(1,0));

        const SafePtr<DGVertex>& ABCD_m = factory.make_child(a,b,c,d,m);
        const SafePtr<DGVertex>& ABCD_mp1 = factory.make_child(a,b,c,d,m+1);
        if (is_simple()) { expr_ = Vector("PA")[dir] * ABCD_m + Vector("WP")[dir] * ABCD_mp1;  nflops_+=3; }

        const F& am1 = a - _1;
        if (exists(am1)) {
          const SafePtr<DGVertex>& Am1BCD_m = factory.make_child(am1,b,c,d,m);
          const SafePtr<DGVertex>& Am1BCD_mp1 = factory.make_child(am1,b,c,d,m+1);
          if (is_simple()) { expr_ += Vector(a)[dir] * Scalar("oo2z") * (Am1BCD_m - Scalar("roz") * Am1BCD_mp1);  nflops_+=5; }
        }
        const F& bm1 = b - _1;
        if (exists(bm1)) {
          const SafePtr<DGVertex>& ABm1CD_m = factory.make_child(a,bm1,c,d,m);
          const SafePtr<DGVertex>& ABm1CD_mp1 = factory.make_child(a,bm1,c,d,m+1);
          if (is_simple()) { expr_ += Vector(b)[dir] * Scalar("oo2z") * (ABm1CD_m - Scalar("roz") * ABm1CD_mp1);  nflops_+=5; }
        }
        const F& cm1 = c - _1;
        if (exists(cm1)) {
          const SafePtr<DGVertex>& ABCm1D_mp1 = factory.make_child(a,b,cm1,d,m+1);
          if (is_simple()) { expr_ += Vector(c)[dir] * Scalar("oo2ze") * ABCm1D_mp1;  nflops_+=3; }
        }
        const F& dm1 = d - _1;
        if (exists(dm1)) {
          const SafePtr<DGVertex>& ABCDm1_mp1 = factory.make_child(a,b,c,dm1,m+1);
          if (is_simple()) { expr_ += Vector(d)[dir] * Scalar("oo2ze") * ABCDm1_mp1;  nflops_+=3; }
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

        const SafePtr<DGVertex>& ABCD_m = factory.make_child(a,b,c,d,m);
        const SafePtr<DGVertex>& ABCD_mp1 = factory.make_child(a,b,c,d,m+1);
        if (is_simple()) { expr_ = Vector("PB")[dir] * ABCD_m + Vector("WP")[dir] * ABCD_mp1;  nflops_+=3; }

        const F& am1 = a - _1;
        if (exists(am1)) {
          const SafePtr<DGVertex>& Am1BCD_m = factory.make_child(am1,b,c,d,m);
          const SafePtr<DGVertex>& Am1BCD_mp1 = factory.make_child(am1,b,c,d,m+1);
          if (is_simple()) { expr_ += Vector(a)[dir] * Scalar("oo2z") * (Am1BCD_m - Scalar("roz") * Am1BCD_mp1);  nflops_+=5; }
        }
        const F& bm1 = b - _1;
        if (exists(bm1)) {
          const SafePtr<DGVertex>& ABm1CD_m = factory.make_child(a,bm1,c,d,m);
          const SafePtr<DGVertex>& ABm1CD_mp1 = factory.make_child(a,bm1,c,d,m+1);
          if (is_simple()) { expr_ += Vector(b)[dir] * Scalar("oo2z") * (ABm1CD_m - Scalar("roz") * ABm1CD_mp1);  nflops_+=5; }
        }
        const F& cm1 = c - _1;
        if (exists(cm1)) {
          const SafePtr<DGVertex>& ABCm1D_mp1 = factory.make_child(a,b,cm1,d,m+1);
          if (is_simple()) { expr_ += Vector(c)[dir] * Scalar("oo2ze") * ABCm1D_mp1;  nflops_+=3; }
        }
        const F& dm1 = d - _1;
        if (exists(dm1)) {
          const SafePtr<DGVertex>& ABCDm1_mp1 = factory.make_child(a,b,c,dm1,m+1);
          if (is_simple()) { expr_ += Vector(d)[dir] * Scalar("oo2ze") * ABCDm1_mp1;  nflops_+=3; }
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

        const SafePtr<DGVertex>& ABCD_m = factory.make_child(a,b,c,d,m);
        const SafePtr<DGVertex>& ABCD_mp1 = factory.make_child(a,b,c,d,m+1);
        if (is_simple()) { expr_ = Vector("QC")[dir] * ABCD_m + Vector("WQ")[dir] * ABCD_mp1;  nflops_+=3; }

        const F& cm1 = c - _1;
        if (exists(cm1)) {
          const SafePtr<DGVertex>& ABCm1D_m = factory.make_child(a,b,cm1,d,m);
          const SafePtr<DGVertex>& ABCm1D_mp1 = factory.make_child(a,b,cm1,d,m+1);
          if (is_simple()) { expr_ += Vector(c)[dir] * Scalar("oo2e") * (ABCm1D_m - Scalar("roe") * ABCm1D_mp1);  nflops_+=5; }
        }
        const F& dm1 = d - _1;
        if (exists(dm1)) {
          const SafePtr<DGVertex>& ABCDm1_m = factory.make_child(a,b,c,dm1,m);
          const SafePtr<DGVertex>& ABCDm1_mp1 = factory.make_child(a,b,c,dm1,m+1);
          if (is_simple()) { expr_ += Vector(d)[dir] * Scalar("oo2e") * (ABCDm1_m - Scalar("roe") * ABCDm1_mp1);  nflops_+=5; }
        }
        const F& am1 = a - _1;
        if (exists(am1)) {
          const SafePtr<DGVertex>& Am1BCD_mp1 = factory.make_child(am1,b,c,d,m+1);
          if (is_simple()) { expr_ += Vector(a)[dir] * Scalar("oo2ze") * Am1BCD_mp1;  nflops_+=3; }
        }
        const F& bm1 = b - _1;
        if (exists(bm1)) {
          const SafePtr<DGVertex>& ABm1CD_mp1 = factory.make_child(a,bm1,c,d,m+1);
          if (is_simple()) { expr_ += Vector(b)[dir] * Scalar("oo2ze") * ABm1CD_mp1;  nflops_+=3; }
        }
        return;
      }
      // Build on D
      if (part == 1 && where == InKet) {
        F a(Tint->bra(0,0));
        F b(Tint->ket(0,0));
        F c(Tint->bra(1,0));
        F d(Tint->ket(1,0) - _1);
        if (!exists(d)) return;

        const SafePtr<DGVertex>& ABCD_m = factory.make_child(a,b,c,d,m);
        const SafePtr<DGVertex>& ABCD_mp1 = factory.make_child(a,b,c,d,m+1);
        if (is_simple()) { expr_ = Vector("QD")[dir] * ABCD_m + Vector("WQ")[dir] * ABCD_mp1;  nflops_+=3; }

        const F& cm1 = c - _1;
        if (exists(cm1)) {
          const SafePtr<DGVertex>& ABCm1D_m = factory.make_child(a,b,cm1,d,m);
          const SafePtr<DGVertex>& ABCm1D_mp1 = factory.make_child(a,b,cm1,d,m+1);
          if (is_simple()) { expr_ += Vector(c)[dir] * Scalar("oo2e") * (ABCm1D_m - Scalar("roe") * ABCm1D_mp1);  nflops_+=5; }
        }
        const F& dm1 = d - _1;
        if (exists(dm1)) {
          const SafePtr<DGVertex>& ABCDm1_m = factory.make_child(a,b,c,dm1,m);
          const SafePtr<DGVertex>& ABCDm1_mp1 = factory.make_child(a,b,c,dm1,m+1);
          if (is_simple()) { expr_ += Vector(d)[dir] * Scalar("oo2e") * (ABCDm1_m - Scalar("roe") * ABCDm1_mp1);  nflops_+=5; }
        }
        const F& am1 = a - _1;
        if (exists(am1)) {
          const SafePtr<DGVertex>& Am1BCD_mp1 = factory.make_child(am1,b,c,d,m+1);
          if (is_simple()) { expr_ += Vector(a)[dir] * Scalar("oo2ze") * Am1BCD_mp1;  nflops_+=3; }
        }
        const F& bm1 = b - _1;
        if (exists(bm1)) {
          const SafePtr<DGVertex>& ABm1CD_mp1 = factory.make_child(a,bm1,c,d,m+1);
          if (is_simple()) { expr_ += Vector(b)[dir] * Scalar("oo2ze") * ABm1CD_mp1;  nflops_+=3; }
        }
        return;
      }
      return;
    }

#if LIBINT_ENABLE_GENERIC_CODE
  template <class F, int part, FunctionPosition where>
    bool
    VRR_11_TwoPRep_11<F,part,where>::has_generic(const SafePtr<CompilationParameters>& cparams) const
    {
      F sh_a(target_->bra(0,0));
      F sh_b(target_->ket(0,0));
      F sh_c(target_->bra(1,0));
      F sh_d(target_->ket(1,0));
      const unsigned int max_opt_am = cparams->max_am_opt();
      if (TrivialBFSet<F>::result)
        return false;
      // generic code works for a0c0 of 0a0c classes where am(a) > 1 and am(c) > 1
      // to generate optimized code for xxxx integral need to generate specialized code for up to (x+x)0(x+x)0 integrals
      if (sh_b.zero() && sh_d.zero() &&
          sh_a.norm() > std::max(2*max_opt_am,1u) && sh_c.norm() > std::max(2*max_opt_am,1u))
        return true;
      if (sh_a.zero() && sh_c.zero() &&
          sh_b.norm() > std::max(2*max_opt_am,1u) && sh_d.norm() > std::max(2*max_opt_am,1u))
        return true;
      return false;
    }

  template <class F, int part, FunctionPosition where>
    std::string
    VRR_11_TwoPRep_11<F,part,where>::generic_header() const
    {
      F sh_a(target_->bra(0,0));
      F sh_b(target_->ket(0,0));
      F sh_c(target_->bra(1,0));
      F sh_d(target_->ket(1,0));
      const bool xsxs = sh_b.zero() && sh_d.zero();
      const bool sxsx = sh_a.zero() && sh_c.zero();
      if (xsxs) return std::string("OSVRR_xs_xs.h");
      if (sxsx) return std::string("OSVRR_sx_sx.h");
      abort(); // unreachable
    }

  template <class F, int part, FunctionPosition where>
    std::string
    VRR_11_TwoPRep_11<F,part,where>::generic_instance(const SafePtr<CodeContext>& context, const SafePtr<CodeSymbols>& args) const
    {
      std::ostringstream oss;
      F sh_a(target_->bra(0,0));
      F sh_b(target_->ket(0,0));
      F sh_c(target_->bra(1,0));
      F sh_d(target_->ket(1,0));
      const bool xsxs = sh_b.zero() && sh_d.zero();
      const bool sxsx = sh_a.zero() && sh_c.zero();

      oss << "using namespace libint2;" << endl;

      if(xsxs) {
        oss << "libint2::OSVRR_xs_xs<" << part << "," << sh_a.norm() << "," << sh_c.norm() << ",";
        oss << ((context->cparams()->max_vector_length() == 1) ? "false" : "true");
        oss << ">::compute(inteval";
      }
      if (sxsx) {
        oss << "libint2::OSVRR_sx_sx<" << part << "," << sh_b.norm() << "," << sh_d.norm() << ",";
        oss << ((context->cparams()->max_vector_length() == 1) ? "false" : "true");
        oss << ">::compute(inteval";
      }

      const unsigned int nargs = args->n();
      for(unsigned int a=0; a<nargs; a++) {
        oss << "," << args->symbol(a);
      }
      oss << ");";

      return oss.str();
    }
#endif // #if !LIBINT_ENABLE_GENERIC_CODE

};

#endif
