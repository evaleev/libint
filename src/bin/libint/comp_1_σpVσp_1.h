//
// Created by Eduard Valeyev on 10/30/23.
//

#ifndef LIBINT_COMP_1_ΣPVΣP_1_H
#define LIBINT_COMP_1_ΣPVΣP_1_H

#include <generic_rr.h>

namespace libint2 {

/**
   * this computes integral of
   * \f$ \sigma \cdot \hat{p} V \sigma \cdot \hat{p} \f$ over CGShell/CGF
   * by rewriting it as a linear combination of integrals over electrostatic
   * potential \f$ V \f$
   * @tparam F basis function type. valid choices are CGShell or CGF
 */
template <typename F>
class CR_1_σpVσp_1 : public GenericRecurrenceRelation< CR_1_σpVσp_1<F>,
                                                    F,
                                                    GenIntegralSet_1_1<F,σpVσpOper,EmptySet> >
{
public:
  typedef CR_1_σpVσp_1<F> ThisType;
  typedef F BasisFunctionType;
  typedef σpVσpOper OperType;
  typedef GenIntegralSet_1_1<F,σpVσpOper,EmptySet> TargetType;
  typedef GenericRecurrenceRelation<ThisType,BasisFunctionType,TargetType> ParentType;
  friend class GenericRecurrenceRelation<ThisType,BasisFunctionType,TargetType>;
  static const unsigned int max_nchildren = 100;  // TODO figure out

  using ParentType::Instance;

  static bool directional() { return false; }

private:
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;
  using ParentType::target_;
  using ParentType::is_simple;

  /// Constructor is private, used by ParentType::Instance that maintains registry of these objects
  CR_1_σpVσp_1(const SafePtr<TargetType>&, unsigned int =0);

  static std::string descr() { return "CR"; }
};

template <typename F>
CR_1_σpVσp_1<F>::CR_1_σpVσp_1(const SafePtr< TargetType >& Tint,
                              unsigned int) : ParentType(Tint,0) {
  assert(Tint->num_func_bra(/* particle */ 0) == 1);
  assert(Tint->num_func_ket(/* particle */ 0) == 1);
  const auto &a = Tint->bra(0, 0);
  const auto &b = Tint->ket(0, 0);
  const auto &oper = Tint->oper();

  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using libint2::algebra::operator*;

  const mType zero_m(0u);

  ChildFactory<ThisType,
               GenIntegralSet_1_1<BasisFunctionType, ElecPotOper, mType>>
      factory(this);

  F Dx_a{a}; Dx_a.deriv().inc(1);
  F Dx_b{b}; Dx_b.deriv().inc(1);
  auto Dx_a_V_Dx_b = factory.make_child(Dx_a, Dx_b, zero_m);
  if (is_simple()) {
    expr_ = Dx_a_V_Dx_b;  // smth like this
    nflops_ += 0;
  }
} // CR_1_σpVσp_1<F>::CR_1_σpVσp_1

}; // namespace libint2

#endif // LIBINT_COMP_1_ΣPVΣP_1_H
