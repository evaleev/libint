
#ifndef _libint2_src_bin_libint_compderivgauss_h_
#define _libint2_src_bin_libint_compderivgauss_h_

#include <generic_rr.h>

using namespace std;

namespace libint2 {

  /** compute relation for derivative 2-e ERI. part+where specify the function
   * to be differentiated.
  */
    template <class IntType, int part, FunctionPosition where>
      class CR_DerivGauss : public GenericRecurrenceRelation< CR_DerivGauss<IntType,part,where>,
                                                              typename IntType::BasisFunctionType,
                                                              IntType >
    {
    public:
      typedef CR_DerivGauss ThisType;
      typedef typename IntType::BasisFunctionType BasisFunctionType;
      typedef IntType TargetType;
      typedef GenericRecurrenceRelation<ThisType,BasisFunctionType,TargetType> ParentType;
      friend class GenericRecurrenceRelation<ThisType,BasisFunctionType,TargetType>;
      static const unsigned int max_nchildren = 2;

      using ParentType::Instance;

      /// always directional! Cartesian derivatives are applied in a particular direction
      static bool directional() { return true; }

    private:
      using ParentType::RecurrenceRelation::expr_;
      using ParentType::RecurrenceRelation::nflops_;
      using ParentType::target_;
      using ParentType::is_simple;

      /// Constructor is private, used by ParentType::Instance that mainains registry of these objects
      CR_DerivGauss(const SafePtr<TargetType>&, unsigned int dir);

      static std::string descr() { return "CR_DerivGauss"; }
      /** Re-Implementation of GenericRecurrenceRelation::generate_label():
          CR Deriv relations are independent of m (it never appears anywhere in equations), hence
          to avoid generating identical code make sure that the (unique) label has m=0. */
      std::string generate_label() const
      {
        typedef typename TargetType::AuxIndexType mType;
        static SafePtr<mType> aux0(new mType(0u));
        ostringstream os;
        os << descr() << "P" << part << to_string(where)
           << genintegralset_label(target_->bra(),target_->ket(),aux0,target_->oper());
        return os.str();
      }
    };

  template <class IntType, int part, FunctionPosition where>
  CR_DerivGauss<IntType,part,where>::CR_DerivGauss(const SafePtr< TargetType >& Tint,
                                                   unsigned int dir) :
    ParentType(Tint,dir)
    {
      using namespace libint2::algebra;
      using namespace libint2::prefactor;
      using namespace libint2::braket;
      typedef BasisFunctionType F;
      const F& _1 = unit<F>(is_simple() ? dir : 0);  // for shell sets increment the first index

      const typename IntType::AuxQuantaType& aux = Tint->aux();
      const typename IntType::OperType& oper = Tint->oper();

      // WARNING assuming one function per position
      { // can't apply to contracted basis functions
        if (where == InBra && Tint->bra(part,0).contracted())
          return;
        if (where == InKet && Tint->ket(part,0).contracted())
          return;
      }

      // the Gaussian must be differentiated in direction dir
      {
        if (where == InBra && Tint->bra(part,0).deriv().d(dir) == 0)
          return;
        if (where == InKet && Tint->ket(part,0).deriv().d(dir) == 0)
          return;
      }

      typedef typename IntType::BraType IBraType;
      typedef typename IntType::KetType IKetType;
      IBraType* bra = new IBraType(Tint->bra());
      IKetType* ket = new IKetType(Tint->ket());

      if (where == InBra) {
        F a(bra->member(part,0));

        // add a+1
        F ap1(bra->member(part,0) + _1);
        ap1.deriv().dec(dir);
        bra->set_member(ap1,part,0);
        auto int_ap1 = this->add_child(IntType::Instance(*bra,*ket,aux,oper));
        bra->set_member(a,part,0);
        if (is_simple()) {
          std::ostringstream oss;
          oss << "two_alpha" << part << "_bra";
          expr_ = Scalar(oss.str()) * int_ap1;  nflops_+=1;
        }

        // See if a-1 exists
        F am1(bra->member(part,0) - _1);
        if (exists(am1)) {
          am1.deriv().dec(dir);
          bra->set_member(am1,part,0);
          auto int_am1 = this->add_child(IntType::Instance(*bra,*ket,aux,oper));
          bra->set_member(a,part,0);
          if (is_simple()) {
            expr_ -= Vector(a)[dir] * int_am1;  nflops_+=2;
          }
        }
        return;
      }

      if (where == InKet) {
        F a(ket->member(part,0));

        // add a+1
        F ap1(ket->member(part,0) + _1);
        ap1.deriv().dec(dir);
        ket->set_member(ap1,part,0);
        auto int_ap1 = this->add_child(IntType::Instance(*bra,*ket,aux,oper));
        ket->set_member(a,part,0);
        if (is_simple()) {
          std::ostringstream oss;
          oss << "two_alpha" << part << "_ket";
          expr_ = Scalar(oss.str()) * int_ap1;  nflops_+=1;
        }

        // See if a-1 exists
        F am1(ket->member(part,0) - _1);
        if (exists(am1)) {
          am1.deriv().dec(dir);
          ket->set_member(am1,part,0);
          auto int_am1 = this->add_child(IntType::Instance(*bra,*ket,aux,oper));
          ket->set_member(a,part,0);
          if (is_simple()) {
            expr_ -= Vector(a)[dir] * int_am1;  nflops_+=2;
          }
        }
        return;
      }

      return;
    }

};

#endif
