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

#ifndef LIBINT_ONEBODY_TASK_TRAITS_H
#define LIBINT_ONEBODY_TASK_TRAITS_H

#include <global_macros.h>
#include <libint2/cgshell_ordering.h>
#include <libint2/config.h>
#include <libint2/shgshell_ordering.h>
#include <master_ints_list.h>

#include <boost/preprocessor/list/for_each.hpp>
#include <cassert>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace libint2 {

/// The operator of the integrals of the 1-body task of operator \p OperType
/// (overlap is computed as the Cartesian multipole of order 0)
template <typename OperType>
using OnebodyTaskOper =
    typename std::conditional<std::is_same<OperType, OverlapOper>::value,
                              CartesianMultipoleOper<3u>, OperType>::type;

/// the defaults of OnebodyTaskTraits
struct OnebodyTaskTraitsBase {
  /// the max order of geometric derivatives the compiler can generate
  static constexpr unsigned int max_deriv_order =
      std::numeric_limits<unsigned int>::max();
  /// true for the operators of the (scalar and X2C) core Hamiltonian, whose
  /// geometric derivatives are generated even if those of the property
  /// operators are disabled (see LIBINT2_DISABLE_ONEBODY_PROPERTY_DERIVS)
  static constexpr bool core_hamiltonian = false;
  /// false if the task is disabled by a configuration option
  static constexpr bool enabled = true;
};

/**
 * Describes how the compiler generates the 1-body task of operator
 * \p OperType (see build_onebody_1b_1k() in build_libint.cc): the static
 * members of OnebodyTaskTraitsBase, and
 * @code
 * using Descr = OnebodyTaskOper<OperType>::Descriptor;
 * /// the operator components of the task with label `label`, in the order of
 * /// the target shell sets (i.e. of Engine::results())
 * static std::vector<Descr> components(const std::string& label);
 * @endcode
 *
 * Specialize it for a new operator rather than testing for the operator type
 * in the code generator.
 */
template <typename OperType>
struct OnebodyTaskTraits;

/// overlap = Cartesian multipole of order 0
template <>
struct OnebodyTaskTraits<OverlapOper> : OnebodyTaskTraitsBase {
  using Descr = CartesianMultipoleOper<3u>::Descriptor;
  static std::vector<Descr> components(const std::string&) { return {Descr()}; }
  static constexpr bool core_hamiltonian = true;
};

template <>
struct OnebodyTaskTraits<KineticOper> : OnebodyTaskTraitsBase {
  using Descr = KineticOper::Descriptor;
  static std::vector<Descr> components(const std::string&) { return {Descr()}; }
  static constexpr bool core_hamiltonian = true;
};

template <>
struct OnebodyTaskTraits<ElecPotOper> : OnebodyTaskTraitsBase {
  using Descr = ElecPotOper::Descriptor;
  static std::vector<Descr> components(const std::string&) { return {Descr()}; }
  static constexpr bool core_hamiltonian = true;
};

/// the Cartesian multipoles of orders 0 through N, where the task label is
/// "Nemultipole"
template <>
struct OnebodyTaskTraits<CartesianMultipoleOper<3u>> : OnebodyTaskTraitsBase {
  using Descr = CartesianMultipoleOper<3u>::Descriptor;
  static std::vector<Descr> components(const std::string& label) {
    // parse the label ... 1emultipole means include multipoles of order 0
    // (overlap) and 1 (dipole)
    const auto key_pos = label.find("emultipole");
    assert(key_pos != std::string::npos);
    unsigned int max_multipole_order = 0;
    std::istringstream iss(label.substr(0, key_pos));
    iss >> max_multipole_order;
    assert(max_multipole_order > 0);

    std::vector<Descr> descrs;
    for (int multipole_order = 0; multipole_order <= max_multipole_order;
         ++multipole_order) {
      // we iterate over them same way as over cartesian Gaussian shells
      int x, y, z;
      LIBINT2_FOR_CART(x, y, z, multipole_order)
      Descr descr;
      descr.inc(0, x);
      descr.inc(1, y);
      descr.inc(2, z);
      descrs.push_back(descr);
      LIBINT2_END_FOR_CART
    }
    return descrs;
  }
};

/// the spherical multipoles of orders 0 through LIBINT_MULTIPOLE_MAX_ORDER
template <>
struct OnebodyTaskTraits<SphericalMultipoleOper> : OnebodyTaskTraitsBase {
  using Descr = SphericalMultipoleOper::Descriptor;
  static std::vector<Descr> components(const std::string&) {
    std::vector<Descr> descrs;
    for (int l = 0; l <= LIBINT_MULTIPOLE_MAX_ORDER; ++l) {
      // we iterate over them using the *standard* solid harmonic Gaussian
      // ordering
      int m;
      LIBINT2_FOR_SOLIDHARM_STANDARD(l, m)
      descrs.emplace_back(l, m);
      LIBINT2_END_FOR_SOLIDHARM
    }
    return descrs;
  }
  /// derivatives of spherical multipole integrals are not implemented
  static constexpr unsigned int max_deriv_order = 0;
};

/// the components of σpVσp: the 4 Pauli components
template <>
struct OnebodyTaskTraits<σpVσpOper> : OnebodyTaskTraitsBase {
  using Descr = σpVσpOper::Descriptor;
  static std::vector<Descr> components(const std::string&) {
    return {Descr(0), Descr(1), Descr(2), Descr(3)};
  }
  static constexpr bool core_hamiltonian = true;
};

/// the components of σpeμσp: 3 dipole directions × 4 Pauli components
template <>
struct OnebodyTaskTraits<σpeμσpOper> : OnebodyTaskTraitsBase {
  using Descr = σpeμσpOper::Descriptor;
  static std::vector<Descr> components(const std::string&) {
    std::vector<Descr> descrs;
    for (int p = 0; p != 12; ++p) descrs.emplace_back(p);
    return descrs;
  }
};

/// the operators with 3 Cartesian components, with descriptors constructed
/// from the axis (x=0, y=1, z=2)
template <typename OperType>
struct AxisVectorOnebodyTaskTraits : OnebodyTaskTraitsBase {
  using Descr = typename OperType::Descriptor;
  static std::vector<Descr> components(const std::string&) {
    return {Descr(0), Descr(1), Descr(2)};
  }
};

/// the magnetic operators (BOOST_PP_MAGNETIC_1_1_OPER_LIST) are generated only
/// if LIBINT2_ENABLE_ONEBODY_MAGNETIC is ON
#ifdef LIBINT_ENABLE_ONEBODY_MAGNETIC
#define LIBINT_ONEBODY_MAGNETIC_ENABLED true
#else
#define LIBINT_ONEBODY_MAGNETIC_ENABLED false
#endif
#define BOOST_PP_MAGNETIC_TASK_TRAITS(r, data, op)                             \
  template <>                                                                  \
  struct OnebodyTaskTraits<op##Oper> : AxisVectorOnebodyTaskTraits<op##Oper> { \
    static constexpr bool enabled = LIBINT_ONEBODY_MAGNETIC_ENABLED;           \
    /* the Engine computes the derivatives of GIAOElecPot w.r.t. the */        \
    /* point charges for the first order only, see */                          \
    /* operator_traits<Operator::giao_nuclear>::max_deriv_order */             \
    static constexpr unsigned int max_deriv_order =                            \
        std::is_same<op##Oper, GIAOElecPotOper>::value                         \
            ? 1                                                                \
            : std::numeric_limits<unsigned int>::max();                        \
  };
BOOST_PP_LIST_FOR_EACH(BOOST_PP_MAGNETIC_TASK_TRAITS, _,
                       BOOST_PP_MAGNETIC_1_1_OPER_LIST)
#undef BOOST_PP_MAGNETIC_TASK_TRAITS
#undef LIBINT_ONEBODY_MAGNETIC_ENABLED

/// @return true if the 1-body task of operator \p OperType with derivative
/// order \p deriv_order is generated
template <typename OperType>
constexpr bool onebody_task_enabled(unsigned int deriv_order) {
  using Traits = OnebodyTaskTraits<OperType>;
  return Traits::enabled && deriv_order <= Traits::max_deriv_order
#ifdef LIBINT_DISABLE_ONEBODY_PROPERTY_DERIVS
         && (deriv_order == 0 || Traits::core_hamiltonian)
#endif
      ;
}

}  // namespace libint2

#endif  // LIBINT_ONEBODY_TASK_TRAITS_H
