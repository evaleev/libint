
#ifndef _libint2_src_bin_libint_cr11g12tig1211_h_
#define _libint2_src_bin_libint_cr11g12tig1211_h_

#include <generic_rr.h>
#include <gaussoper.h>

using namespace std;

namespace libint2 {

  /** Compute relation for 2-e integrals of the G12_Ti_G12 operators.
  */
  template <class BFSet>
    class CR_11_G12TiG12_11 : public GenericRecurrenceRelation< CR_11_G12TiG12_11<BFSet>,
                                                             BFSet,
                                                             GenIntegralSet_11_11<BFSet,G12TiG12,mType> >
    {
    public:
      typedef CR_11_G12TiG12_11 ThisType;
      typedef BFSet BasisFunctionType;
      typedef GenIntegralSet_11_11<BFSet,G12TiG12,mType> TargetType;
      typedef GenericRecurrenceRelation<ThisType,BFSet,TargetType> ParentType;
      friend class GenericRecurrenceRelation<ThisType,BFSet,TargetType>;
      static const unsigned int max_nchildren = 1;

      using ParentType::Instance;
    private:
      using ParentType::RecurrenceRelation::expr_;
      using ParentType::RecurrenceRelation::nflops_;
      using ParentType::target_;
      using ParentType::is_simple;

      /// Constructor is private, used by ParentType::Instance that mainains registry of these objects
      CR_11_G12TiG12_11(const SafePtr<TargetType>&, unsigned int dir);
      /// This relation is not directional
      static bool directional() { return false; }
      static std::string descr() { return "CR"; }

#if LIBINT_ENABLE_GENERIC_CODE
    /// Implementation of RecurrenceRelation::has_generic()
    bool has_generic(const SafePtr<CompilationParameters>& cparams) const;
    /// Implementation of RecurrenceRelation::generic_header()
    std::string generic_header() const { return "GenericScale.h"; }
    /// Implementation of RecurrenceRelation::generic_instance()
    std::string generic_instance(const SafePtr<CodeContext>& context, const SafePtr<CodeSymbols>& args) const;
#endif
    };

  template <class F>
    CR_11_G12TiG12_11<F>::CR_11_G12TiG12_11(const SafePtr<TargetType>& Tint,
                                            unsigned int dir) :
    ParentType(Tint,dir)
    {
      if (dir != 0)
        return;
      using namespace libint2::algebra;
      using namespace libint2::prefactor;
      using namespace libint2::braket;
      // kinetic energy of which electron?
      const int i = target_->oper()->descr().K();
      const R12kG12 G2(2);

      F a(Tint->bra(0,0));
      F b(Tint->ket(0,0));
      F c(Tint->bra(1,0));
      F d(Tint->ket(1,0));

      if (target_->oper()->descr().contracted())
        return;

      // [G12,[T1,G12]] = [G12,[T2,G12]]
      {
        typedef GenIntegralSet_11_11<BasisFunctionType,R12kG12,mType> ChildType;
        ChildFactory<ThisType,ChildType> factory(this);
        auto ab_G2_cd = factory.make_child(a,b,c,d,0u,G2);
        if (is_simple()) {
          expr_ = Scalar("R12_2_G12_scale_to_G12T1G12") * ab_G2_cd;
          nflops_ += 1;
        }
      }

    }

#if LIBINT_ENABLE_GENERIC_CODE
  template <class F>
  bool
  CR_11_G12TiG12_11<F>::has_generic(const SafePtr<CompilationParameters>& cparams) const
  {
    F sh_a(target_->bra(0,0));
    F sh_b(target_->ket(0,0));
    F sh_c(target_->bra(1,0));
    F sh_d(target_->ket(1,0));
    const unsigned int max_opt_am = cparams->max_am_opt();
    // to generate optimized code for xxxx integral need to generate specialized code for up to (x+x)0(x+x)0 integrals
    if (!TrivialBFSet<F>::result &&
        (sh_a.norm() > max_opt_am ||
         sh_b.norm() > max_opt_am ||
         sh_c.norm() > max_opt_am ||
         sh_d.norm() > max_opt_am
        )
       )
      return true;
    return false;
  }

  template <class F>
  std::string
  CR_11_G12TiG12_11<F>::generic_instance(const SafePtr<CodeContext>& context, const SafePtr<CodeSymbols>& args) const {
      std::ostringstream oss;

      const bool vec = (context->cparams()->max_vector_length() != 1);
      if (vec)
        oss << "_libint2_static_api_scale_vec_short_(";
      else
        oss << "_libint2_static_api_scale_short_(";

      const unsigned int nargs = args->n();
      for(unsigned int a=0; a<nargs; a++) {
        oss << args->symbol(a) << ",";
      }

      oss << target_->size() << "*" << (vec ? "inteval->veclen" : "1")
          << ",inteval->R12_2_G12_scale_to_G12T1G12" << (vec ? ",inteval->veclen" : "[0]") << ");";

      return oss.str();
  }
#endif

};

#endif
