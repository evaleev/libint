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

#ifndef _libint2_src_bin_libint_oper_h_
#define _libint2_src_bin_libint_oper_h_

#include <contractable.h>
#include <global_macros.h>
#include <hashable.h>
#include <iter.h>
#include <multipole.h>
#include <util.h>
#include <vectorn.h>

#include <boost/preprocessor/list/for_each.hpp>
#include <string>

namespace libint2 {

/** Permutational symmetries: antisymmetric(anti), symmetric(symm), nonsymmetric
   (nonsymm), some more complicated symmetry (nonstd) */
struct PermutationalSymmetry {
  typedef enum { anti = -1, symm = 1, nonsymm = 0, nonstd = -2 } type;
};

/** OperatorProperties describes various properties of an operator or operator
   set
    @tparam NP number of particles
    @tparam multi true if multiplicative
    @tparam psymmetry symmetry with respect to permutation of bra and ket
    @tparam origin_dependent true if operator is origin-dependent
  */
template <unsigned int NP, bool multi, PermutationalSymmetry::type psymmetry,
          bool origin_dependent = false>
class OperatorProperties {
 public:
  static constexpr auto np = NP;
  static constexpr auto multiplicative = multi;
  static constexpr auto psymm = psymmetry;
  static constexpr auto odep = origin_dependent;
};

/** OperSet is the base class for all (sets of) operators.
   OperSet's must be constructable using
   std::shared_ptr<OperSet> or std::shared_ptr<ConstructablePolymorphically>.
*/
class OperSet : public ConstructablePolymorphically {
 public:
  typedef DummyIterator iter_type;
  virtual ~OperSet(){};

  /// Returns full human-readable description of the operator
  virtual std::string description() const = 0;
  /// Returns short label for the operator
  virtual std::string label() const = 0;
  /** Returns 1, 0, or -1, if each operator in the set is symmetric,
    nonsymmetric, or antisymmetric with respect to permutation of particles i
    and j */
  virtual int psymm(int i, int j) const = 0;
  /** Returns 1, 0, or -1, if each operator in the set is Hermitian,
    non-Hermitian, or anti-Hermitian w.r.t. particle p */
  virtual int hermitian(int p) const = 0;

  /// is operator origin-dependent?
  virtual bool origin_dependent() const = 0;

  /// Number of operators in the set
  virtual unsigned int num_oper() const = 0;
};

/** Oper is OperSet characterized by properties Props
 */
template <class Props>
class Oper : public OperSet {
 public:
  typedef Props Properties;
  virtual ~Oper() {}

  /// Implementation of OperSet::psymm()
  int psymm(int i, int j) const override;
  /// Implementation of OperSet::hermitian()
  int hermitian(int p) const override;
  /// Implementation of OperSet::origin_dependent()
  bool origin_dependent() const override { return Props::odep; }

  bool operator==(const Oper&) const;

 protected:
  /// The only declared constructor is only useable by derived classes
  Oper() {}

 private:
  /// Must be overloaded by derived classes if Props::psymm ==
  /// PermutationalSymmetry::nonstd
  virtual int nonstd_psymm(int i, int j) const {
    throw ProgrammingError("nonstd_psymm is not overloaded");
  }
  /// Must be overloaded by derived classes if Props::multi != true
  virtual int nonstd_hermitian(int p) const {
    throw ProgrammingError("nonstd_hermitian is not overloaded");
  }
};

template <class Props>
int Oper<Props>::psymm(int i, int j) const {
  if (i < 0 || i >= static_cast<int>(Props::np))
    throw std::runtime_error(
        "Oper<Props>::psymm(i,j) -- index i out of bounds");
  if (j < 0 || j >= static_cast<int>(Props::np))
    throw std::runtime_error(
        "Oper<Props>::psymm(i,j) -- index j out of bounds");
  if (i == j) return 1;

  switch (Props::psymm) {
    case PermutationalSymmetry::anti:
      return -1;
    case PermutationalSymmetry::symm:
      return 1;
    case PermutationalSymmetry::nonsymm:
      return 0;
    case PermutationalSymmetry::nonstd:
      return nonstd_psymm(i, j);
    default:
      abort();
  }
}

template <class Props>
int Oper<Props>::hermitian(int p) const {
  if (Props::multiplicative)
    return +1;
  else
    return nonstd_hermitian(p);
}

template <class Props>
bool Oper<Props>::operator==(const Oper& a) const {
  return true;
}

//////////////////////////////

/** GenOper is a single operator described by descriptor Descr
 */
template <class Descr>
class GenOper : public Oper<typename Descr::Properties>,
                public Hashable<unsigned, ComputeKey> {
 public:
  typedef Descr Descriptor;
  typedef typename Descr::Properties Properties;
  typedef Oper<Properties> parent_type;
  /// GenOper is not a set
  typedef GenOper iter_type;

  unsigned int num_oper() const override { return 1; }
  /// Implementation of Hashable::key()
  unsigned int key() const override { return descr_.key(); }
  /// Range of key is [0,Descr::max_key)
  static const unsigned int max_key = Descr::max_key;
  /// Implementation of OperSet::description()
  std::string description() const override { return descr_.description(); }
  /// Implementation of OperSet::label()
  std::string label() const override { return descr_.label(); }
  /// Return the descriptor object
  Descr& descr() { return descr_; }
  /// Return the descriptor object
  const Descr& descr() const { return descr_; }

  GenOper(Descr descr = Descr()) : descr_(descr) {}
  GenOper(const std::shared_ptr<GenOper>& o) : descr_(o->descr_) {}
  GenOper(const std::shared_ptr<OperSet>& o)
      : descr_(require_dynamic_cast<GenOper, OperSet>(o)->descr_) {}
  GenOper(const std::shared_ptr<ConstructablePolymorphically>& o)
      : descr_(require_dynamic_cast<GenOper, ConstructablePolymorphically>(o)
                   ->descr_) {}
  explicit GenOper(const ConstructablePolymorphically& o)
      : descr_(require_dynamic_cast<GenOper, ConstructablePolymorphically>(&o)
                   ->descr_) {}
  virtual ~GenOper() {}

 private:
  Descr descr_;

  /// Implementation of Oper::nonstd_psymm()
  int nonstd_psymm(int i, int j) const override {
    // TODO: figure out how to call this only if Desc::Properties::psymm ==
    // PermutationalSymmetry::nonstd
    if (Descr::Properties::psymm == PermutationalSymmetry::nonstd)
      return descr_.psymm(i, j);
    else
      throw ProgrammingError(
          "GenOper::nonstd_psymm -- descriptor is not nonstd");
  }

  /// Implementation of Oper::nonstd_hermitian()
  int nonstd_hermitian(int i) const override {
    // TODO: figure out how to call this only if Desc::Properties::psymm ==
    // PermutationalSymmetry::nonstd
    if (!Descr::Properties::multiplicative)
      return descr_.hermitian(i);
    else
      throw ProgrammingError(
          "GenOper::nonstd_hermitian -- this operator is multiplicative");
  }
};

//////////////////////////////

typedef OperatorProperties<1, false, PermutationalSymmetry::nonsymm>
    Nonmultiplicative1Body_Props;
typedef OperatorProperties<1, true, PermutationalSymmetry::nonsymm>
    Multiplicative1Body_Props;
typedef OperatorProperties<1, true, PermutationalSymmetry::nonsymm, true>
    MultiplicativeODep1Body_Props;
typedef OperatorProperties<2, true, PermutationalSymmetry::symm>
    MultiplicativeSymm2Body_Props;
typedef OperatorProperties<2, true, PermutationalSymmetry::nonsymm>
    MultiplicativeNonsymm2Body_Props;
typedef OperatorProperties<2, false, PermutationalSymmetry::symm>
    NonmultiplicativeSymm2Body_Props;
typedef OperatorProperties<2, false, PermutationalSymmetry::nonsymm>
    NonmultiplicativeNonsymm2Body_Props;

/** GenMultSymmOper is a generic multiplicative symmetric N-body operator
 */
template <unsigned int N>
struct GenMultSymmOper_Descr : public Contractable<GenMultSymmOper_Descr<N>> {
  typedef OperatorProperties<N, true, PermutationalSymmetry::symm> Properties;
  static const unsigned int max_key = 1;
  unsigned int key() const { return 0; }
  std::string description() const {
    return "generic multiplicative symmetric operator";
  }
  std::string label() const { return "GenMultSymmOper"; }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { abort(); }
};
typedef GenOper<GenMultSymmOper_Descr<2>> GenMultSymm2BodyOper;

#define BOOST_PP_DECLARE_HERMITIAN_ONEBODY_DESCRIPTOR(r, propprefix, opname) \
  struct opname##_Descr : public Contractable<opname##_Descr> {              \
    typedef propprefix##1Body_Props Properties;                              \
    static const unsigned int max_key = 1;                                   \
    unsigned int key() const { return 0; }                                   \
    std::string description() const { return #opname; }                      \
    std::string label() const { return #opname; }                            \
    int psymm(int i, int j) const { abort(); }                               \
    int hermitian(int i) const { return +1; }                                \
  };                                                                         \
  typedef GenOper<opname##_Descr> opname##Oper;

#define BOOST_PP_HERMITIAN_ONEBODY_OPER_LIST (Kinetic, BOOST_PP_NIL)
BOOST_PP_LIST_FOR_EACH(BOOST_PP_DECLARE_HERMITIAN_ONEBODY_DESCRIPTOR,
                       Nonmultiplicative, BOOST_PP_HERMITIAN_ONEBODY_OPER_LIST)
#undef BOOST_PP_HERMITIAN_ONEBODY_OPER_LIST
#define BOOST_PP_HERMITIAN_ONEBODY_OPER_LIST (Overlap, BOOST_PP_NIL)
BOOST_PP_LIST_FOR_EACH(BOOST_PP_DECLARE_HERMITIAN_ONEBODY_DESCRIPTOR,
                       Multiplicative, BOOST_PP_HERMITIAN_ONEBODY_OPER_LIST)
#undef BOOST_PP_HERMITIAN_ONEBODY_OPER_LIST
#define BOOST_PP_HERMITIAN_ONEBODY_OPER_LIST (ElecPot, BOOST_PP_NIL)
BOOST_PP_LIST_FOR_EACH(BOOST_PP_DECLARE_HERMITIAN_ONEBODY_DESCRIPTOR,
                       MultiplicativeODep, BOOST_PP_HERMITIAN_ONEBODY_OPER_LIST)
#undef BOOST_PP_HERMITIAN_ONEBODY_OPER_LIST

// N.B. σpVσp has 4 (Pauli) components, so need to customize its descriptor ...
// may need to generalize this if need similar operators elsewhere
// #define BOOST_PP_HERMITIAN_ONEBODY_OPER_LIST (σpVσp, BOOST_PP_NIL)
// BOOST_PP_LIST_FOR_EACH ( BOOST_PP_DECLARE_HERMITIAN_ONEBODY_DESCRIPTOR,
// MultiplicativeODep, BOOST_PP_HERMITIAN_ONEBODY_OPER_LIST) #undef
// BOOST_PP_HERMITIAN_ONEBODY_OPER_LIST

struct σpVσp_Descr : public Contractable<σpVσp_Descr> {
  typedef MultiplicativeODep1Body_Props Properties;

  σpVσp_Descr() : quaternion_index_(0) {}
  σpVσp_Descr(int quaternion_index) : quaternion_index_(quaternion_index) {
    assert(quaternion_index <= 3);
  }

  static const unsigned int max_key = 4;
  unsigned int key() const { return quaternion_index(); }
  std::string description() const {
    std::string descr("opVop[");
    if (quaternion_index() == 0)
      descr += "0";
    else if (quaternion_index() == 1)
      descr += "X";
    else if (quaternion_index() == 2)
      descr += "Y";
    else if (quaternion_index() == 3)
      descr += "Z";
    else
      abort();
    return descr + "]";
  }
  std::string label() const { return description(); }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { return +1; }

  int quaternion_index() const { return quaternion_index_; }

 private:
  const int quaternion_index_ = -1;
};
typedef GenOper<σpVσp_Descr> σpVσpOper;

/** opemuop: (μ σ·p | r | σ·p ν), one-body σ·p-on-both-sides analog of dipole.
 *  σ_a σ_b = δ_ab + iε_abc σ_c folds the 9 raw σ_a∂_a r_k σ_b∂_b dyadics per
 *  dipole direction k down to 4 Pauli-quaternion components (trace + 3
 *  antisym), mirroring σpVσp's fold of σ·p V σ·p. 12 outputs total = 3 dipole
 *  directions × 4 Pauli components, indexed composite_index = 4·k + q.
 */
struct σpeμσp_Descr : public Contractable<σpeμσp_Descr> {
  typedef MultiplicativeODep1Body_Props Properties;

  σpeμσp_Descr() : composite_index_(0) {}
  σpeμσp_Descr(int composite_index) : composite_index_(composite_index) {
    assert(composite_index >= 0 && composite_index < 12);
  }

  static const unsigned int max_key = 12;
  unsigned int key() const { return composite_index(); }
  std::string description() const {
    static const char* dipole_lbl[] = {"x", "y", "z"};
    static const char* pauli_lbl[] = {"0", "X", "Y", "Z"};
    const auto ci = composite_index();
    if (ci < 0 || ci >= 12) abort();
    return std::string("opemuop[") + dipole_lbl[ci / 4] + "," +
           pauli_lbl[ci % 4] + "]";
  }
  std::string label() const { return description(); }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { return +1; }

  int composite_index() const { return composite_index_; }
  /// dipole direction ∈ {0=x, 1=y, 2=z}
  int dipole_dir() const { return composite_index_ / 4; }
  /// Pauli quaternion ∈ {0=trace, 1=σ_x, 2=σ_y, 3=σ_z}
  int quaternion_index() const { return composite_index_ % 4; }

 private:
  const int composite_index_ = -1;
};
typedef GenOper<σpeμσp_Descr> σpeμσpOper;

/// cartesian multipole operator in \c NDIM dimensions
/// \f$ \hat{O}(\vec{k}) \equiv \vec{r}^{\cdot \vec{k}} = r_1^{k_1} r_2^{k_2}
/// \cdots \f$ \internal OriginDerivative<NDIM> is used to store cartesian
/// multipole quanta, not the derivative quanta
template <unsigned int NDIM>
struct CartesianMultipole_Descr
    : public Contractable<CartesianMultipole_Descr<NDIM>>,
      public CartesianMultipoleQuanta<NDIM> {
  typedef MultiplicativeODep1Body_Props Properties;
  using CartesianMultipoleQuanta<NDIM>::max_key;

  CartesianMultipole_Descr() {}
  CartesianMultipole_Descr(unsigned int k) {
    assert(NDIM == 1u);
    this->inc(0, k);
  }
  std::string description() const {
    std::string descr("CartesianMultipole[");
    std::ostringstream oss;
    for (unsigned i = 0; i != NDIM; ++i) {
      oss << (*this)[i];
      if (i + 1 != NDIM) oss << ",";
    }
    return descr + oss.str() + "]";
  }
  std::string label() const { return description(); }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { return +1; }
};
template <unsigned int NDIM>
using CartesianMultipoleOper = GenOper<CartesianMultipole_Descr<NDIM>>;

/** Represents quantum numbers of \em real spherical multipole operator
 * defined in Eqs. 5 and 6 of J.M. Pérez-Jordá and W. Yang, J Chem Phys 104,
 * 8003 (1996). ( \f$ m \geq 0 \f$ corresponds to moments \f$ \mathcal{N}^+ \f$
 * , \f$ m < 0 \f$ corresponds to \f$ \mathcal{N}^- \f$ )
 */
struct SphericalMultipole_Descr : public Contractable<SphericalMultipole_Descr>,
                                  public SphericalMultipoleQuanta {
  typedef MultiplicativeODep1Body_Props Properties;
  using SphericalMultipoleQuanta::max_key;
  using SphericalMultipoleQuanta::Sign;

  /// Default ctor makes a 0th-order multipole
  SphericalMultipole_Descr() : SphericalMultipole_Descr(0, 0) {}
  /// constructs \f$ \mathcal{N}^{+}_{l,m} \f$ if \f$ m \geq 0 \f$, otherwise
  /// constructs \f$ \mathcal{N}^{-}_{l,m} \f$
  SphericalMultipole_Descr(int l, int m) : SphericalMultipoleQuanta(l, m) {}
  SphericalMultipole_Descr(int l, int m, Sign sign)
      : SphericalMultipoleQuanta(l, m, sign) {}
  SphericalMultipole_Descr(const SphericalMultipoleQuanta& quanta)
      : SphericalMultipoleQuanta(quanta) {}

  std::string description() const {
    std::string descr =
        std::string("SphericalMultipole[") + std::to_string(this->l()) + "," +
        std::to_string((this->sign() == Sign::plus ? 1 : -1) * this->m()) + "]";
    return descr;
  }
  std::string label() const { return description(); }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { return +1; }
};
using SphericalMultipoleOper = GenOper<SphericalMultipole_Descr>;

/** TwoPRep is the two-body repulsion operator.
 */
struct TwoPRep_Descr : public Contractable<TwoPRep_Descr> {
  typedef MultiplicativeSymm2Body_Props Properties;
  static const unsigned int max_key = 1;
  unsigned int key() const { return 0; }
  std::string description() const { return "1/r_{12}"; }
  std::string label() const { return "TwoPRep"; }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { return +1; }
};
typedef GenOper<TwoPRep_Descr> TwoPRep;

/** Coulombσpσp is the two-body repulsion operator.
 */
struct Coulombσpσp_Descr : public Contractable<Coulombσpσp_Descr> {
  typedef MultiplicativeSymm2Body_Props Properties;

  Coulombσpσp_Descr() : quaternion_index_(0) {}
  Coulombσpσp_Descr(int quaternion_index)
      : quaternion_index_(quaternion_index) {
    assert(quaternion_index <= 3);
  }

  static const unsigned int max_key = 4;
  unsigned int key() const { return quaternion_index(); }
  std::string description() const {
    std::string descr("coulomb_opop[");
    if (quaternion_index() == 0)
      descr += "0";
    else if (quaternion_index() == 1)
      descr += "X";
    else if (quaternion_index() == 2)
      descr += "Y";
    else if (quaternion_index() == 3)
      descr += "Z";
    else
      abort();
    return descr + "]";
  }
  std::string label() const { return description(); }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { return +1; }

  int quaternion_index() const { return quaternion_index_; }

 private:
  const int quaternion_index_ = -1;
};
typedef GenOper<Coulombσpσp_Descr> CoulombσpσpOper;

/** Coulombσp is the 3-center single-σ·p operator:
 *  (P | 1/r_{12} | μ, σ·p ν), with the fitting function P a spectator and a
 *  single σ·p acting on the 2nd ket AO function ν. Unlike Coulombσpσp (two σ·p
 *  on the ket pair, folded via the Dirac identity into 4 quaternion
 *  components), a single σ·p produces no σ·σ folding: the result is the bare
 *  Cartesian gradient of ν, i.e. the SO(3) vector irrep with 3 components
 *  (x=0, y=1, z=2). Each component is exactly one first-derivative ERI child;
 *  the imaginary unit and Pauli matrix of σ·p are applied by the caller.
 */
struct Coulombσp_Descr : public Contractable<Coulombσp_Descr> {
  typedef MultiplicativeSymm2Body_Props Properties;

  Coulombσp_Descr() : cartesian_index_(0) {}
  Coulombσp_Descr(int cartesian_index) : cartesian_index_(cartesian_index) {
    assert(cartesian_index >= 0 && cartesian_index <= 2);
  }

  static const unsigned int max_key = 3;
  unsigned int key() const { return cartesian_index(); }
  std::string description() const {
    std::string descr("coulomb_op[");
    if (cartesian_index() == 0)
      descr += "X";
    else if (cartesian_index() == 1)
      descr += "Y";
    else if (cartesian_index() == 2)
      descr += "Z";
    else
      abort();
    return descr + "]";
  }
  std::string label() const { return description(); }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { return +1; }

  int cartesian_index() const { return cartesian_index_; }

 private:
  const int cartesian_index_ = -1;
};
typedef GenOper<Coulombσp_Descr> CoulombσpOper;

struct σpσpCoulombσpσp_Descr : public Contractable<σpσpCoulombσpσp_Descr> {
  typedef MultiplicativeSymm2Body_Props Properties;

  σpσpCoulombσpσp_Descr() : quaternion_index_(0) {}
  σpσpCoulombσpσp_Descr(int quaternion_index)
      : quaternion_index_(quaternion_index) {
    assert(quaternion_index >= 0 && quaternion_index <= 15);
  }

  /// 16 components from tensor product of two independent spin spaces:
  /// index = 4 * bra_spin_index + ket_spin_index
  /// where spin indices are: 0=S (scalar), 1=X, 2=Y, 3=Z (cross product)
  static const unsigned int max_key = 16;
  unsigned int key() const { return quaternion_index(); }
  std::string description() const {
    // clang-format off
    // Option A (tensor product order): index = 4 * bra_spin + ket_spin
    static const char* labels[] = {
        "SS", "SX", "SY", "SZ",
        "XS", "XX", "XY", "XZ",
        "YS", "YX", "YY", "YZ",
        "ZS", "ZX", "ZY", "ZZ"
    };
    // clang-format on
    const auto qi = quaternion_index();
    if (qi > 15) abort();
    return std::string("opop_coulomb_opop[") + labels[qi] + "]";
  }
  std::string label() const { return description(); }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { return +1; }

  int quaternion_index() const { return quaternion_index_; }

 private:
  const int quaternion_index_ = -1;
};
typedef GenOper<σpσpCoulombσpσp_Descr> σpσpCoulombσpσpOper;

/** σpCoulombσp: (μ σ·p ν | 1/r_{12} | κ σ·p λ).
 *  Gaunt LS "bilinear" operator with one σ·p on each side.
 *  Outputs the SO(3) irreducible decomposition of the 3×3 gradient-gradient
 *  tensor T_{ab} = ∂_a ∂_b (μν|κλ): 1 scalar trace + 3 antisymmetric
 *  (curl-curl) + 5 symmetric-traceless = 9 components. Same storage as the
 *  raw dyadic, but indexed by physics-meaningful irreps so consumers do not
 *  need to hand-build trace/antisym/sym-TL combinations at every contraction
 *  site.
 */
struct σpCoulombσp_Descr : public Contractable<σpCoulombσp_Descr> {
  typedef MultiplicativeSymm2Body_Props Properties;

  /// SO(3) irreducible components of the rank-2 Cartesian tensor T_{ab}.
  enum Component : int {
    Scalar = 0,      ///< T_xx + T_yy + T_zz (trace, ∇·∇)
    AntisymX = 1,    ///< T_yz − T_zy = (∇×∇)_x
    AntisymY = 2,    ///< T_zx − T_xz = (∇×∇)_y
    AntisymZ = 3,    ///< T_xy − T_yx = (∇×∇)_z
    SymTLDiagA = 4,  ///< T_xx − T_yy
    SymTLDiagB = 5,  ///< 2·T_zz − T_xx − T_yy
    SymTLOffXY = 6,  ///< T_xy + T_yx
    SymTLOffXZ = 7,  ///< T_xz + T_zx
    SymTLOffYZ = 8,  ///< T_yz + T_zy
  };

  σpCoulombσp_Descr() : component_index_(0) {}
  σpCoulombσp_Descr(int component_index) : component_index_(component_index) {
    assert(component_index >= 0 && component_index <= 8);
  }

  static const unsigned int max_key = 9;
  unsigned int key() const { return component_index(); }
  std::string description() const {
    static const char* labels[] = {
        "scalar",       "antisym_x",    "antisym_y",
        "antisym_z",    "symtl_diag_a", "symtl_diag_b",
        "symtl_off_xy", "symtl_off_xz", "symtl_off_yz",
    };
    const auto ci = component_index();
    if (ci > 8) abort();
    return std::string("op_coulomb_op[") + labels[ci] + "]";
  }
  std::string label() const { return description(); }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { return +1; }

  int component_index() const { return component_index_; }

  bool is_scalar() const { return component_index_ == Scalar; }
  bool is_antisym() const {
    return component_index_ >= AntisymX && component_index_ <= AntisymZ;
  }
  bool is_sym_tl() const {
    return component_index_ >= SymTLDiagA && component_index_ <= SymTLOffYZ;
  }
  /// for antisym components: 0=x, 1=y, 2=z (only valid if is_antisym())
  int antisym_cart() const { return component_index_ - AntisymX; }

 private:
  const int component_index_ = -1;
};
typedef GenOper<σpCoulombσp_Descr> σpCoulombσpOper;

/** GTG_1d is the two-body 1-dimensional Gaussian geminal
 */
struct GTG_1d_Descr : public Contractable<GTG_1d_Descr> {
  typedef MultiplicativeSymm2Body_Props Properties;
  static const unsigned int max_key = 1;
  unsigned int key() const { return 0; }
  std::string description() const { return "GTG_1d"; }
  std::string label() const { return "GTG_1d"; }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { return +1; }
};
typedef GenOper<GTG_1d_Descr> GTG_1d;

/** R12_k_G12 is a two-body operator of form r_{12}^k * exp(-\gamma * r_{12}),
    where k is an integer and \gamma is a positive real number.
*/
class R12_k_G12_Descr : public Contractable<R12_k_G12_Descr> {
 public:
  typedef MultiplicativeSymm2Body_Props Properties;
  /// K can range from -1 to 4
  R12_k_G12_Descr(int K) : K_(K) { assert(K >= -1 && K <= 4); }
  R12_k_G12_Descr(const R12_k_G12_Descr& a)
      : Contractable<R12_k_G12_Descr>(a), K_(a.K_) {}
  static const unsigned int max_key = 5;
  unsigned int key() const { return K_ + 1; }
  std::string description() const { return label_(K_, this->contracted()); }
  std::string label() const { return symbol_(K_, this->contracted()); }
  int K() const { return K_; }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { abort(); }

 private:
  R12_k_G12_Descr();
  static std::string label_(int K, bool contracted);
  static std::string symbol_(int K, bool contracted);
  int K_;
};
typedef GenOper<R12_k_G12_Descr> R12kG12;

/** R12k_R12l_G12 is a two-body operator of form ( r_{12x}^kx * r_{12y}^ky *
   r_{12z}^kz ) * (r_{12x}^lx * r_{12y}^ly * r_{12z}^lz ) * G12 The following
   restrictions are imposed: 0 <= kx+ky+kz <= 4, 0 <= lx+ly+lz <= 4
  */
class R12k_R12l_G12_Descr : public Contractable<R12k_R12l_G12_Descr> {
 public:
  typedef MultiplicativeSymm2Body_Props Properties;
  static const int kmax = 4;
  R12k_R12l_G12_Descr(const IntVec3& K, const IntVec3& L) : K_(K), L_(L) {}
  R12k_R12l_G12_Descr(const R12k_R12l_G12_Descr& a)
      : Contractable<R12k_R12l_G12_Descr>(a), K_(a.K_), L_(a.L_) {}
  const IntVec3& K() const { return K_; }
  const IntVec3& L() const { return L_; }
  static const unsigned int max_key = kmax * kmax * kmax * kmax * kmax * kmax;
  unsigned int key() const;
  std::string description() const { return label_(K_, L_, this->contracted()); }
  std::string label() const { return symbol_(K_, L_, this->contracted()); }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { abort(); }

 private:
  R12k_R12l_G12_Descr();
  static std::string label_(const IntVec3& K, const IntVec3& L,
                            bool contracted);
  static std::string symbol_(const IntVec3& K, const IntVec3& L,
                             bool contracted);
  IntVec3 K_;
  IntVec3 L_;
};
typedef GenOper<R12k_R12l_G12_Descr> R12kR12lG12;

/** Ti_G12 is a two-body operator of form [T_i, G12],
    where i is particle index (0 or 1) and G12 is a Gaussian Geminal.
*/
class Ti_G12_Descr : public Contractable<Ti_G12_Descr> {
 public:
  typedef NonmultiplicativeNonsymm2Body_Props Properties;
  /// K can range from 0 to 1
  static const unsigned int max_key = 2;
  Ti_G12_Descr(int K) : K_(K) { assert(K >= 0 && K <= 1); }
  Ti_G12_Descr(const Ti_G12_Descr& a)
      : Contractable<Ti_G12_Descr>(a), K_(a.K_) {}
  unsigned int key() const { return K_; }
  std::string description() const { return label_(K_, this->contracted()); }
  std::string label() const { return symbol_(K_, this->contracted()); }
  int K() const { return K_; }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const {
    if (i != K_)
      return +1;
    else
      return -1;
  }

 private:
  Ti_G12_Descr();
  static std::string label_(int K, bool contracted);
  static std::string symbol_(int K, bool contracted);
  int K_;
};
typedef GenOper<Ti_G12_Descr> TiG12;

/** G12_Ti_G12 is a two-body operator of form [G12, [T_i, G12]] = (Nabla_i G12)
   \dot (Nabla_i G12) where i is particle index (0 or 1) and G12 is a Gaussian
   Geminal. It is a *multiplicative* operator!
*/
class G12_Ti_G12_Descr : public Contractable<G12_Ti_G12_Descr> {
 public:
  typedef MultiplicativeSymm2Body_Props Properties;
  /// K can range from 0 to 1
  static const unsigned int max_key = 2;
  G12_Ti_G12_Descr(int K) : K_(K) { assert(K >= 0 && K <= 1); }
  G12_Ti_G12_Descr(const G12_Ti_G12_Descr& a)
      : Contractable<G12_Ti_G12_Descr>(a), K_(a.K_) {}
  unsigned int key() const { return K_; }
  std::string description() const { return label_(K_, this->contracted()); }
  std::string label() const { return symbol_(K_, this->contracted()); }
  int K() const { return K_; }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { return +1; }

 private:
  G12_Ti_G12_Descr();
  static std::string label_(int K, bool contracted);
  static std::string symbol_(int K, bool contracted);
  int K_;
};
typedef GenOper<G12_Ti_G12_Descr> G12TiG12;

/** r_1.r_1 x g12 is a result of differentiation of exp( - a r_1^2 - a r_2^2 - c
 * r_{12}^2) geminal .
 */
struct R1dotR1_G12_Descr : public Contractable<R1dotR1_G12_Descr> {
  typedef MultiplicativeNonsymm2Body_Props Properties;
  static const unsigned int max_key = 1;
  unsigned int key() const { return 0; }
  std::string description() const { return "r_1.r_1 x G12"; }
  std::string label() const { return "R1dotR1_G12"; }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { abort(); }
};
typedef GenOper<R1dotR1_G12_Descr> R1dotR1_G12;

/** r_2.r_2 x g12 is a result of differentiation of exp( - a r_1^2 - a r_2^2 - c
 * r_{12}^2) geminal .
 */
struct R2dotR2_G12_Descr : public Contractable<R2dotR2_G12_Descr> {
  typedef MultiplicativeNonsymm2Body_Props Properties;
  static const unsigned int max_key = 1;
  unsigned int key() const { return 0; }
  std::string description() const { return "r_2.r_2 x G12"; }
  std::string label() const { return "R2dotR2_G12"; }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { abort(); }
};
typedef GenOper<R2dotR2_G12_Descr> R2dotR2_G12;

/** r_1.r_2 x g12 is a result of differentiation of exp( - a r_1^2 - a r_2^2 - c
 * r_{12}^2) geminal .
 */
struct R1dotR2_G12_Descr : public Contractable<R1dotR2_G12_Descr> {
  typedef MultiplicativeSymm2Body_Props Properties;
  static const unsigned int max_key = 1;
  unsigned int key() const { return 0; }
  std::string description() const { return "r_1.r_2 x G12"; }
  std::string label() const { return "R1dotR2_G12"; }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { abort(); }
};
typedef GenOper<R1dotR2_G12_Descr> R1dotR2_G12;

/** \f$ (\nabla_I \cdot g_{12}') (g_{12}' \cdot \nabla_K) \f$ is a component of
   \f$ [g_{12}, [\nabla_I^4, g_{12}]] \f$ integral where I = 1 or 2.
*/
struct DivG12prime_xTx_Descr : public Contractable<DivG12prime_xTx_Descr> {
  typedef NonmultiplicativeNonsymm2Body_Props Properties;
  static const unsigned int max_key = 2;
  DivG12prime_xTx_Descr(int I) : I_(I) { assert(I >= 0 && I <= 1); }
  DivG12prime_xTx_Descr(const DivG12prime_xTx_Descr& a)
      : Contractable<DivG12prime_xTx_Descr>(a), I_(a.I_) {}
  unsigned int key() const { return I_; }
  std::string description() const { return label_(I_); }
  std::string label() const { return symbol_(I_); }
  int I() const { return I_; }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const {
    if (i != I_)
      return +1;
    else
      return -1;
  }

 private:
  DivG12prime_xTx_Descr();
  static std::string label_(int I);
  static std::string symbol_(int I);
  int I_;
};
typedef GenOper<DivG12prime_xTx_Descr> DivG12prime_xTx;

};  // namespace libint2

#endif
