
#ifndef _libint2_src_bin_libint_compderivgauss_h_
#define _libint2_src_bin_libint_compderivgauss_h_

#include <generic_rr.h>

using namespace std;

namespace libint2 {

  /** VRR Recurrence Relation for 2-e ERI. part specifies for which particle
  the angular momentum is raised. where specifies whether angular momentum is decreased in bra or ket.
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
    private:
      using ParentType::RecurrenceRelation::expr_;
      using ParentType::RecurrenceRelation::nflops_;
      using ParentType::target_;
      using ParentType::is_simple;

      /// Constructor is private, used by ParentType::Instance that mainains registry of these objects
      CR_DerivGauss(const SafePtr<TargetType>&, unsigned int dir);
      /// Default directionality
      static bool directional() { return ParentType::default_directional(); }

      static std::string descr() { return "CR_DerivGauss"; }
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
