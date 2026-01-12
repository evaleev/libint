/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint compiler.
 *
 *  Libint compiler is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint compiler is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint compiler.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <comp_xyz.h>
#include <dg.h>
#include <graph_registry.h>
#include <integral_1_1.h>
#include <intset_to_ints.h>
#include <prefactors.h>
#include <rr.h>
#include <singl_stack.h>
#include <strategy.h>
#include <uncontract.h>

#include <boost/preprocessor/list/for_each.hpp>

using namespace std;

namespace libint2 {

template <>
void CR_XYZ_1_1<CGShell, OverlapOper, EmptySet>::compute(const CGShell& a,
                                                         const CGShell& b,
                                                         const OperType&) {
  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using namespace libint2::braket;

  // why do I need this? when this function is in header clang-602.0.49 on OS X
  // works without this
  using libint2::algebra::operator*;

  ChildFactory<ThisType, GenIntegralSet_1_1<CGShell1d<CartesianAxis_X>,
                                            OverlapOper, EmptySet>>
      factory_X(this);
  auto ab_X = factory_X.make_child(a.norm() + a.deriv().norm(),
                                   b.norm() + b.deriv().norm());

  ChildFactory<ThisType, GenIntegralSet_1_1<CGShell1d<CartesianAxis_Y>,
                                            OverlapOper, EmptySet>>
      factory_Y(this);
  auto ab_Y = factory_Y.make_child(a.norm() + a.deriv().norm(),
                                   b.norm() + b.deriv().norm());

  ChildFactory<ThisType, GenIntegralSet_1_1<CGShell1d<CartesianAxis_Z>,
                                            OverlapOper, EmptySet>>
      factory_Z(this);
  auto ab_Z = factory_Z.make_child(a.norm() + a.deriv().norm(),
                                   b.norm() + b.deriv().norm());
}

template <>
void CR_XYZ_1_1<CGShell1d<CartesianAxis_X>, OverlapOper, EmptySet>::compute(
    const CGShell1d<CartesianAxis_X>& a, const CGShell1d<CartesianAxis_X>& b,
    const OperType&) {
  this->add_child(prefactor::Scalar("_0_Overlap_0_x"));
}
template <>
void CR_XYZ_1_1<CGShell1d<CartesianAxis_Y>, OverlapOper, EmptySet>::compute(
    const CGShell1d<CartesianAxis_Y>& a, const CGShell1d<CartesianAxis_Y>& b,
    const OperType&) {
  this->add_child(prefactor::Scalar("_0_Overlap_0_y"));
}
template <>
void CR_XYZ_1_1<CGShell1d<CartesianAxis_Z>, OverlapOper, EmptySet>::compute(
    const CGShell1d<CartesianAxis_Z>& a, const CGShell1d<CartesianAxis_Z>& b,
    const OperType&) {
  this->add_child(prefactor::Scalar("_0_Overlap_0_z"));
}

template <>
void CR_XYZ_1_1<CGF, OverlapOper, EmptySet>::compute(const CGF& a, const CGF& b,
                                                     const OperType&) {
  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using namespace libint2::braket;

  // why do I need this? when this function is in header clang-602.0.49 on OS X
  // works without this
  using libint2::algebra::operator*;

  ChildFactory<ThisType, GenIntegralSet_1_1<CGF1d<CartesianAxis_X>, OverlapOper,
                                            EmptySet>>
      factory_X(this);
  auto ab_X = factory_X.make_child(a[CartesianAxis_X], b[CartesianAxis_X]);

  ChildFactory<ThisType, GenIntegralSet_1_1<CGF1d<CartesianAxis_Y>, OverlapOper,
                                            EmptySet>>
      factory_Y(this);
  auto ab_Y = factory_Y.make_child(a[CartesianAxis_Y], b[CartesianAxis_Y]);

  ChildFactory<ThisType, GenIntegralSet_1_1<CGF1d<CartesianAxis_Z>, OverlapOper,
                                            EmptySet>>
      factory_Z(this);
  auto ab_Z = factory_Z.make_child(a[CartesianAxis_Z], b[CartesianAxis_Z]);

  if (is_simple()) {
    expr_ = ab_X * ab_Y * ab_Z;
    nflops_ += 2;
  }

}  // end of Overlap compute

template <>
void CR_XYZ_1_1<CGShell, KineticOper, EmptySet>::compute(const CGShell& a,
                                                         const CGShell& b,
                                                         const OperType&) {
  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using namespace libint2::braket;

  ChildFactory<ThisType, GenIntegralSet_1_1<CGShell1d<CartesianAxis_X>,
                                            OverlapOper, EmptySet>>
      factory_X(this);
  auto ab_X = factory_X.make_child(a.norm() + a.deriv().norm() + 1,
                                   b.norm() + b.deriv().norm() + 1);

  ChildFactory<ThisType, GenIntegralSet_1_1<CGShell1d<CartesianAxis_Y>,
                                            OverlapOper, EmptySet>>
      factory_Y(this);
  auto ab_Y = factory_Y.make_child(a.norm() + a.deriv().norm() + 1,
                                   b.norm() + b.deriv().norm() + 1);

  ChildFactory<ThisType, GenIntegralSet_1_1<CGShell1d<CartesianAxis_Z>,
                                            OverlapOper, EmptySet>>
      factory_Z(this);
  auto ab_Z = factory_Z.make_child(a.norm() + a.deriv().norm() + 1,
                                   b.norm() + b.deriv().norm() + 1);
}

template <>
void CR_XYZ_1_1<CGF, KineticOper, EmptySet>::compute(const CGF& a, const CGF& b,
                                                     const OperType&) {
  using namespace libint2::algebra;
  // why do I need this? when this function is in header clang-602.0.49 on OS X
  // works without this
  using libint2::algebra::operator*;

  // this is what I want this code to look like
  // this->wedge(dot(Nabla(a), Nabla(b)), Scalar(-0.25), EmptySet(),
  // OverlapOper());

  ChildFactory<ThisType, GenIntegralSet_1_1<CGF1d<CartesianAxis_X>, OverlapOper,
                                            EmptySet>>
      factory_X(this);
  auto s_X = factory_X.make_child(a[CartesianAxis_X], b[CartesianAxis_X]);

  ChildFactory<ThisType, GenIntegralSet_1_1<CGF1d<CartesianAxis_Y>, OverlapOper,
                                            EmptySet>>
      factory_Y(this);
  auto s_Y = factory_Y.make_child(a[CartesianAxis_Y], b[CartesianAxis_Y]);

  ChildFactory<ThisType, GenIntegralSet_1_1<CGF1d<CartesianAxis_Z>, OverlapOper,
                                            EmptySet>>
      factory_Z(this);
  auto s_Z = factory_Z.make_child(a[CartesianAxis_Z], b[CartesianAxis_Z]);

  ChildFactory<ThisType, GenIntegralSet_1_1<CGF1d<CartesianAxis_X>, KineticOper,
                                            EmptySet>>
      tfactory_X(this);
  auto t_X = tfactory_X.make_child(a[CartesianAxis_X], b[CartesianAxis_X]);

  ChildFactory<ThisType, GenIntegralSet_1_1<CGF1d<CartesianAxis_Y>, KineticOper,
                                            EmptySet>>
      tfactory_Y(this);
  auto t_Y = tfactory_Y.make_child(a[CartesianAxis_Y], b[CartesianAxis_Y]);

  ChildFactory<ThisType, GenIntegralSet_1_1<CGF1d<CartesianAxis_Z>, KineticOper,
                                            EmptySet>>
      tfactory_Z(this);
  auto t_Z = tfactory_Z.make_child(a[CartesianAxis_Z], b[CartesianAxis_Z]);

  expr_ = t_X * s_Y * s_Z;
  expr_ += s_X * t_Y * s_Z;
  expr_ += s_X * s_Y * t_Z;
  nflops_ += 8;
}  // end of Kinetic 3-d compute

#define BOOST_PP_CR_XYZ_1_1_1D_KINETIC_COMPUTE(r, data, elem)                  \
  template <>                                                                  \
  void                                                                         \
  CR_XYZ_1_1<CGF1d<CartesianAxis_##elem>, KineticOper, EmptySet>::compute(     \
      const BasisFunctionType& a, const BasisFunctionType& b,                  \
      const OperType&) {                                                       \
    using namespace libint2::algebra;                                          \
    using namespace libint2::prefactor;                                        \
    using libint2::algebra::operator*;                                         \
                                                                               \
    ChildFactory<ThisType,                                                     \
                 GenIntegralSet_1_1<BasisFunctionType, OverlapOper, EmptySet>> \
        factory(this);                                                         \
                                                                               \
    const BasisFunctionType& _1 = unit<BasisFunctionType>(0);                  \
    auto ap1 = a + _1;                                                         \
    auto bp1 = b + _1;                                                         \
    auto ap1bp1 = factory.make_child(ap1, bp1);                                \
    expr_ = Scalar(0.5) * Scalar("two_alpha0_bra") *                           \
            Scalar("two_alpha0_ket") * ap1bp1;                                 \
    nflops_ += 3;                                                              \
    auto am1 = a - _1;                                                         \
    auto bm1 = b - _1;                                                         \
    if (exists(am1)) {                                                         \
      auto am1bp1 = factory.make_child(am1, bp1);                              \
      expr_ -= Scalar(0.5 * a[0]) * Scalar("two_alpha0_ket") * am1bp1;         \
      nflops_ += 3;                                                            \
    }                                                                          \
    if (exists(bm1)) {                                                         \
      auto ap1bm1 = factory.make_child(ap1, bm1);                              \
      expr_ -= Scalar(0.5 * b[0]) * Scalar("two_alpha0_bra") * ap1bm1;         \
      nflops_ += 3;                                                            \
    }                                                                          \
    if (exists(am1) && exists(bm1)) {                                          \
      auto am1bm1 = factory.make_child(am1, bm1);                              \
      expr_ += Scalar(0.5 * a[0] * b[0]) * am1bm1;                             \
      nflops_ += 2;                                                            \
    }                                                                          \
  }

#define BOOST_PP_XYZ_LIST (X, (Y, (Z, BOOST_PP_NIL)))
BOOST_PP_LIST_FOR_EACH(BOOST_PP_CR_XYZ_1_1_1D_KINETIC_COMPUTE, _,
                       BOOST_PP_XYZ_LIST)
#undef BOOST_PP_XYZ_LIST

template <>
void CR_XYZ_1_1<CGShell, CartesianMultipoleOper<3u>, EmptySet>::compute(
    const CGShell& a, const CGShell& b, const OperType& oper) {
  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using namespace libint2::braket;

  // why do I need this? when this function is in header clang-602.0.49 on OS X
  // works without this
  using libint2::algebra::operator*;

  ChildFactory<ThisType, GenIntegralSet_1_1<CGShell1d<CartesianAxis_X>,
                                            OverlapOper, EmptySet>>
      factory_X(this);
  auto ab_X =
      factory_X.make_child(a.norm() + a.deriv().norm(),
                           b.norm() + b.deriv().norm() + oper.descr().norm());

  ChildFactory<ThisType, GenIntegralSet_1_1<CGShell1d<CartesianAxis_Y>,
                                            OverlapOper, EmptySet>>
      factory_Y(this);
  auto ab_Y =
      factory_Y.make_child(a.norm() + a.deriv().norm(),
                           b.norm() + b.deriv().norm() + oper.descr().norm());

  ChildFactory<ThisType, GenIntegralSet_1_1<CGShell1d<CartesianAxis_Z>,
                                            OverlapOper, EmptySet>>
      factory_Z(this);
  auto ab_Z =
      factory_Z.make_child(a.norm() + a.deriv().norm(),
                           b.norm() + b.deriv().norm() + oper.descr().norm());
}

template <>
void CR_XYZ_1_1<CGF, CartesianMultipoleOper<3u>, EmptySet>::compute(
    const CGF& a, const CGF& b, const OperType& oper) {
  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using namespace libint2::braket;

  // why do I need this? when this function is in header clang-602.0.49 on OS X
  // works without this
  using libint2::algebra::operator*;

  using Descr1d = CartesianMultipoleOper<1u>::Descriptor;

  ChildFactory<ThisType,
               GenIntegralSet_1_1<CGF1d<CartesianAxis_X>,
                                  CartesianMultipoleOper<1u>, EmptySet>>
      factory_X(this);
  auto ab_X =
      factory_X.make_child(a[CartesianAxis_X], b[CartesianAxis_X], EmptySet(),
                           Descr1d(oper.descr()[CartesianAxis_X]));

  ChildFactory<ThisType,
               GenIntegralSet_1_1<CGF1d<CartesianAxis_Y>,
                                  CartesianMultipoleOper<1u>, EmptySet>>
      factory_Y(this);
  auto ab_Y =
      factory_Y.make_child(a[CartesianAxis_Y], b[CartesianAxis_Y], EmptySet(),
                           Descr1d(oper.descr()[CartesianAxis_Y]));

  ChildFactory<ThisType,
               GenIntegralSet_1_1<CGF1d<CartesianAxis_Z>,
                                  CartesianMultipoleOper<1u>, EmptySet>>
      factory_Z(this);
  auto ab_Z =
      factory_Z.make_child(a[CartesianAxis_Z], b[CartesianAxis_Z], EmptySet(),
                           Descr1d(oper.descr()[CartesianAxis_Z]));

  expr_ = ab_X * ab_Y * ab_Z;
  nflops_ += 2;
}

#define BOOST_PP_CR_XYZ_1_1_1D_CMULTIPOLE_COMPUTE(r, data, elem)              \
  template <>                                                                 \
  void CR_XYZ_1_1<CGF1d<CartesianAxis_##elem>, CartesianMultipoleOper<1u>,    \
                  EmptySet>::compute(const BasisFunctionType& a,              \
                                     const BasisFunctionType& b,              \
                                     const OperType& oper) {                  \
    using namespace libint2::algebra;                                         \
    using namespace libint2::prefactor;                                       \
    using libint2::algebra::operator*;                                        \
                                                                              \
    const BasisFunctionType& _1 = unit<BasisFunctionType>(0);                 \
    auto bp1 = b + _1;                                                        \
    auto descr_m1 = oper.descr();                                             \
    if (descr_m1[0] > 0) {                                                    \
      ChildFactory<ThisType,                                                  \
                   GenIntegralSet_1_1<BasisFunctionType, OperType, EmptySet>> \
          factory(this);                                                      \
      descr_m1.dec(0);                                                        \
      auto oper_m1 = OperType(descr_m1);                                      \
      expr_ = factory.make_child(a, bp1, EmptySet(), oper_m1) +               \
              Vector("BO")[CartesianAxis_##elem] *                            \
                  factory.make_child(a, b, EmptySet(), oper_m1);              \
    } else {                                                                  \
      ChildFactory<ThisType, GenIntegralSet_1_1<BasisFunctionType,            \
                                                OverlapOper, EmptySet>>       \
          sfactory(this);                                                     \
      expr_ = Scalar(0u) + sfactory.make_child(a, b);                         \
    }                                                                         \
    nflops_ += 2;                                                             \
  }

#define BOOST_PP_XYZ_LIST (X, (Y, (Z, BOOST_PP_NIL)))
BOOST_PP_LIST_FOR_EACH(BOOST_PP_CR_XYZ_1_1_1D_CMULTIPOLE_COMPUTE, _,
                       BOOST_PP_XYZ_LIST)
#undef BOOST_PP_XYZ_LIST

}  // namespace libint2
