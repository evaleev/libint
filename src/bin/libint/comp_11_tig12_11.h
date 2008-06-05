
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <dgvertex.h>
#include <rr.h>
#include <integral.h>
#include <tig12_11_11.h>
#include <algebra.h>
#include <flop.h>
#include <prefactors.h>
#include <context.h>
#include <default_params.h>

#ifndef _libint2_src_bin_libint_cr11tig1211_h_
#define _libint2_src_bin_libint_cr11tig1211_h_

using namespace std;


namespace libint2 {

  /** Compute relation for 2-e integrals of the Ti_G12 operators.
  */
  template <template <typename,typename,typename> class I, class BFSet>
  class CR_11_TiG12_11 : public RecurrenceRelation
    {

  public:
    typedef RecurrenceRelation ParentType;
    typedef BFSet BasisFunctionType;
    typedef CR_11_TiG12_11<I,BFSet> ThisType;
    typedef GenIntegralSet_11_11<BFSet,TiG12,mType> TargetType;
    typedef GenIntegralSet_11_11<BFSet,R12kG12,mType> ChildType;
    /// The type of expressions in which RecurrenceRelations result.
    typedef RecurrenceRelation::ExprType ExprType;

    /** Use Instance() to obtain an instance of RR. This function is provided to avoid
        issues with getting a SafePtr from constructor (as needed for registry to work).
        second argument is dummy to make sure that the interface is consistent among all RRs.
    */
    static SafePtr<ThisType> Instance(const SafePtr<TargetType>&, unsigned int);
    ~CR_11_TiG12_11() {}

    /// Implementation of RecurrenceRelation::num_children()
    const unsigned int num_children() const { return children_.size(); };
    /// Implementation of RecurrenceRelation::rr_target()
    SafePtr<DGVertex> rr_target() const { return static_pointer_cast<DGVertex,TargetType>(target_); }
    /// Implementation of RecurrenceRelation::rr_child()
    SafePtr<DGVertex> rr_child(unsigned int i) const { return dynamic_pointer_cast<DGVertex,ChildType>(children_.at(i)); }
    /// Implementation of RecurrenceRelation::is_simple()
    bool is_simple() const {
      return TrivialBFSet<BFSet>::result;
    }
    
    const std::string cpp_function_name() {}
    const std::string cpp_source_name() {}
    const std::string cpp_header_name() {}
    std::ostream& cpp_source(std::ostream&) {}

  private:
    CR_11_TiG12_11(const SafePtr<TargetType>&);

    SafePtr<TargetType> target_;
    static const unsigned int max_nchildren_ = 18;
    std::vector< SafePtr<ChildType> > children_;
    const SafePtr<ChildType>& make_child(const BFSet& A, const BFSet& B, const BFSet& C, const BFSet& D, int K) {
      typedef typename ChildType::OperType OperType;
      // [Ti,G12] is reduced to R12^k * G12
      const OperType oper(K);
      const SafePtr<ChildType>& i = ChildType::Instance(A,B,C,D,mType(0),oper);
      children_.push_back(i);
      return *(children_.end()-1);
    }

    std::string generate_label() const
    {
      ostringstream os;
      os << "CR ( " << rr_target()->label() << " )";
      return os.str();
    }
  };

  template <template <typename,typename,typename> class I, class F>
    SafePtr< CR_11_TiG12_11<I,F> >
    CR_11_TiG12_11<I,F>::Instance(const SafePtr<TargetType>& Tint, unsigned int xyz)
    {
      SafePtr<ThisType> this_ptr(new ThisType(Tint));
      // Do post-construction duties
      if (this_ptr->num_children() != 0) {
        this_ptr->register_with_rrstack<ThisType>();
      }
      return this_ptr;
    }

  template <template <typename,typename,typename> class I, class F>
    CR_11_TiG12_11<I,F>::CR_11_TiG12_11(const SafePtr<TargetType>& Tint) :
    target_(Tint)
    {
      using namespace libint2::algebra;
      using namespace libint2::prefactor;
      children_.reserve(max_nchildren_);
      // kinetic energy of which electron?
      const int i = target_->oper()->descr().K();

      // TODO rederive compute relationship as [Ti,G12] = -1/2 [\nabla,[\nabla,G12]] - [\nabla,G12] \cdot \nabla
      throw std::logic_error("CR_11_TiG12_11 has not yet been re-implemented");
      
      // [T1,G12]
      if (i == 0) {
        
      }
      // [T1,G12]
      if (i == 1) {
        
      }
      
    }
  
};

#endif
