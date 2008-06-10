
#ifndef _libint2_src_bin_libint_cr11tig1211_h_
#define _libint2_src_bin_libint_cr11tig1211_h_

#include <generic_rr.h>
#include <tig12_11_11.h>

using namespace std;

namespace libint2 {
    
  /** Compute relation for 2-e integrals of the Ti_G12 operators.
  */
  template <class BFSet>
    class CR_11_TiG12_11 : public GenericRecurrenceRelation< CR_11_TiG12_11<BFSet>,
                                                             BFSet,
                                                             GenIntegralSet_11_11<BFSet,TiG12,mType>,
                                                             GenIntegralSet_11_11<BFSet,R12kG12,mType> >
    {
    public:
      typedef CR_11_TiG12_11 ThisType;
      typedef BFSet BasisFunctionType;
      typedef GenIntegralSet_11_11<BFSet,TiG12,mType> TargetType;
      typedef GenIntegralSet_11_11<BFSet,R12kG12,mType> ChildType;
      typedef GenericRecurrenceRelation<ThisType,BFSet,TargetType,ChildType> ParentType;
      friend class GenericRecurrenceRelation<ThisType,BFSet,TargetType,ChildType>;
      static const unsigned int max_nchildren = 19;

      using ParentType::Instance;
    private:
      using RecurrenceRelation::expr_;
      using RecurrenceRelation::nflops_;
      using ParentType::target_;
      using ParentType::is_simple;

      /// Constructor is private, used by ParentType::Instance that mainains registry of these objects
      CR_11_TiG12_11(const SafePtr<TargetType>&, unsigned int dir);
      /// helper function that makes writing the contructor easier
      const SafePtr<ChildType>& make_child(const BFSet& A, const BFSet& B, const BFSet& C, const BFSet& D, int K) {
        typedef typename ChildType::OperType OperType;
        // [Ti,G12] is reduced to R12^k * G12
        const OperType oper(K);
        const SafePtr<ChildType>& i = ChildType::Instance(A,B,C,D,mType(0),oper);
        return add_child(i);
      }
      static std::string descr() { return "CR"; }

    };

  template <class F>
    CR_11_TiG12_11<F>::CR_11_TiG12_11(const SafePtr<TargetType>& Tint,
                                      unsigned int dir) :
    ParentType(Tint,dir)
    {
      using namespace libint2::algebra;
      using namespace libint2::prefactor;
      // kinetic energy of which electron?
      const int i = target_->oper()->descr().K();

      // TODO rederive compute relationship as [Ti,G12] = -1/2 [\nabla,[\nabla,G12]] - [\nabla,G12] \cdot \nabla
      throw std::logic_error("CR_11_TiG12_11 has not yet been re-implemented");
      
      // [T1,G12]
      if (i == 0) {
        
      }
      // [T2,G12]
      if (i == 1) {
        
      }
      
    }
  
};

#endif
