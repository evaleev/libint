
#ifndef _libint2_src_bin_libint_cr11divg12primextx11_h_
#define _libint2_src_bin_libint_cr11divg12primextx11_h_

#include <generic_rr.h>
#include <integral_11_11.h>
#include <gaussoper.h>

using namespace std;

namespace libint2 {
    
  /** Compute relation for 2-e integrals of the DivG12prime_xTx operators.
  */
  template <class BFSet>
    class CR_11_DivG12prime_xTx_11 : public GenericRecurrenceRelation< CR_11_DivG12prime_xTx_11<BFSet>,
                                                             BFSet,
                                                             GenIntegralSet_11_11<BFSet,DivG12prime_xTx,mType>,
                                                             GenIntegralSet_11_11<BFSet,R12kG12,mType> >
    {
    public:
      typedef CR_11_DivG12prime_xTx_11 ThisType;
      typedef BFSet BasisFunctionType;
      typedef GenIntegralSet_11_11<BFSet,DivG12prime_xTx,mType> TargetType;
      typedef GenIntegralSet_11_11<BFSet,R12kG12,mType> ChildType;
      typedef GenericRecurrenceRelation<ThisType,BFSet,TargetType,ChildType> ParentType;
      friend class GenericRecurrenceRelation<ThisType,BFSet,TargetType,ChildType>;
      static const unsigned int max_nchildren = 100;

      using ParentType::Instance;
    private:
      using RecurrenceRelation::expr_;
      using RecurrenceRelation::nflops_;
      using ParentType::target_;
      using ParentType::is_simple;

      /// Constructor is private, used by ParentType::Instance that mainains registry of these objects
      CR_11_DivG12prime_xTx_11(const SafePtr<TargetType>&, unsigned int dir);
      static std::string descr() { return "CR"; }

    };

  template <class F>
    CR_11_DivG12prime_xTx_11<F>::CR_11_DivG12prime_xTx_11(const SafePtr<TargetType>& Tint,
                                      unsigned int dir) :
    ParentType(Tint,dir)
    {
      using namespace libint2::algebra;
      using namespace libint2::prefactor;
      using namespace libint2::braket;
      
      F a(Tint->bra(0,0));
      F b(Tint->ket(0,0));
      F c(Tint->bra(1,0));
      F d(Tint->ket(1,0));

      // |Nabla1.G12' ac ) |Nabla1.G12' bd )
      {
        ParentType::wedge( R12vec_dot_Nabla1(_pbra(a,c)) , R12vec_dot_Nabla1(_pket(b,d)), mType(0u), R12kG12(0));
        if (is_simple()) expr_ *= Scalar(4.0) * Scalar("gamma_bra") * Scalar("gamma_ket");
      }

    }
  
};


#endif
