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

  σpVσp_Descr() : pauli_index_(0) {}
  σpVσp_Descr(int pauli_index) : pauli_index_(pauli_index) {
    assert(pauli_index <= 3);
  }

  static const unsigned int max_key = 4;
  unsigned int key() const { return pauli_index(); }
  std::string description() const {
    std::string descr("opVop[");
    if (pauli_index() == 0)
      descr += "0";
    else if (pauli_index() == 1)
      descr += "Z";
    else if (pauli_index() == 2)
      descr += "X";
    else if (pauli_index() == 3)
      descr += "Y";
    else
      abort();
    return descr + "]";
  }
  std::string label() const { return description(); }
  int psymm(int i, int j) const { abort(); }
  int hermitian(int i) const { return +1; }

  int pauli_index() const { return pauli_index_; }

 private:
  const int pauli_index_ = -1;
};
typedef GenOper<σpVσp_Descr> σpVσpOper;

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
