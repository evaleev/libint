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

  // can express integrals of σpVσp in terms of derivative integrals of V for primitive Gaussians only
  if (a.contracted() || b.contracted())
    return;

  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using libint2::algebra::operator*;

  const mType zero_m(0u);

  ChildFactory<ThisType,
               GenIntegralSet_1_1<BasisFunctionType, ElecPotOper, mType>>
      factory(this);

  constexpr auto x = 0;
  constexpr auto y = 1;
  constexpr auto z = 2;

  // (a|W0|b) = (d a/dAx | V | d b/dBx) + (d a/dAy | V | d b/dBy) + (d a/dAz | V | d b/dBz)
  switch (oper->descr().pauli_index()) {
  case 0: {
    F Dx_a{a};
    Dx_a.deriv().inc(x);
    F Dx_b{b};
    Dx_b.deriv().inc(x);
    auto Dx_a_V_Dx_b = factory.make_child(Dx_a, Dx_b, zero_m);
    F Dy_a{a};
    Dy_a.deriv().inc(y);
    F Dy_b{b};
    Dy_b.deriv().inc(y);
    auto Dy_a_V_Dy_b = factory.make_child(Dy_a, Dy_b, zero_m);
    F Dz_a{a};
    Dz_a.deriv().inc(z);
    F Dz_b{b};
    Dz_b.deriv().inc(z);
    auto Dz_a_V_Dz_b = factory.make_child(Dz_a, Dz_b, zero_m);
    if (is_simple()) {
      expr_ = Scalar(-1) *
              (Dx_a_V_Dx_b + Dy_a_V_Dy_b + Dz_a_V_Dz_b);
      nflops_ += 3;
    }
  } break;
  // N.B.: Wz is needed after W0 because it's needed for 2-component methods. Wx and Wy are not.
  // (a|Wz|b) = (d a/dAx | V | d b/dBy) - (d a/dAy | V | d b/dBx)
  case 1: {
    F Dx_a{a};
    Dx_a.deriv().inc(x);
    F Dy_b{b};
    Dy_b.deriv().inc(y);
    auto Dx_a_V_Dy_b = factory.make_child(Dx_a, Dy_b, zero_m);
    F Dy_a{a};
    Dy_a.deriv().inc(y);
    F Dx_b{b};
    Dx_b.deriv().inc(x);
    auto Dy_a_V_Dx_b = factory.make_child(Dy_a, Dx_b, zero_m);
    if (is_simple()) {
      expr_ = Scalar(-1) *
              (Dx_a_V_Dy_b - Dy_a_V_Dx_b);
      nflops_ += 2;
    }
  } break;
  // (a|Wx|b) = (d a/dAy | V | d b/dBz) - (d a/dAz | V | d b/dBy)
  case 2: {
    F Dy_a{a};
    Dy_a.deriv().inc(y);
    F Dz_b{b};
    Dz_b.deriv().inc(z);
    auto Dy_a_V_Dz_b = factory.make_child(Dy_a, Dz_b, zero_m);
    F Dz_a{a};
    Dz_a.deriv().inc(z);
    F Dy_b{b};
    Dy_b.deriv().inc(y);
    auto Dz_a_V_Dy_b = factory.make_child(Dz_a, Dy_b, zero_m);
    if (is_simple()) {
      expr_ = Scalar(-1) *
              (Dy_a_V_Dz_b - Dz_a_V_Dy_b);
      nflops_ += 2;
    }
  } break;
  // (a|Wy|b) = (d a/dAz | V | d b/dBx) - (d a/dAx | V | d b/dBz)
  case 3: {
    F Dz_a{a};
    Dz_a.deriv().inc(z);
    F Dx_b{b};
    Dx_b.deriv().inc(x);
    auto Dz_a_V_Dx_b = factory.make_child(Dz_a, Dx_b, zero_m);
    F Dx_a{a};
    Dx_a.deriv().inc(x);
    F Dz_b{b};
    Dz_b.deriv().inc(z);
    auto Dx_a_V_Dz_b = factory.make_child(Dx_a, Dz_b, zero_m);
    if (is_simple()) {
      expr_ = Scalar(-1) *
              (Dz_a_V_Dx_b - Dx_a_V_Dz_b);
      nflops_ += 2;
    }
  } break;
  default:
    throw std::runtime_error("CR_1_σpVσp_1: invalid Pauli index");
  }

} // CR_1_σpVσp_1<F>::CR_1_σpVσp_1

}; // namespace libint2

#endif // LIBINT_COMP_1_ΣPVΣP_1_H
