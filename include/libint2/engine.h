/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint library.
 *
 *  Libint library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_lib_libint_engine_h_
#define _libint2_src_lib_libint_engine_h_

#ifndef LIBINT2_DOES_NOT_INLINE_ENGINE
#define __libint2_engine_inline inline
#else
#define __libint2_engine_inline
#endif

#include <libint2/util/cxxstd.h>
#if LIBINT2_CPLUSPLUS_STD < 2011
#error "libint2/engine.h requires C++11 support"
#endif

#include <libint2/boys_fwd.h>
#include <libint2/braket.h>
#include <libint2/cartesian.h>
#include <libint2/cgshellinfo.h>
#include <libint2/cxxapi.h>
#include <libint2/shell.h>
#include <libint2/solidharmonics.h>
#include <libint2/util/any.h>
#include <libint2/util/array_adaptor.h>
#include <libint2/util/compressed_pair.h>
#include <libint2/util/intpart_iter.h>
#include <libint2/util/timer.h>

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

// the engine will be profiled by default if library was configured with
// --enable-profile
#ifdef LIBINT2_PROFILE
#define LIBINT2_ENGINE_TIMERS
// uncomment if want to profile each integral class
#define LIBINT2_ENGINE_PROFILE_CLASS
#endif
// uncomment if want to profile the engine even if library was configured
// without --enable-profile
// #  define LIBINT2_ENGINE_TIMERS

namespace libint2 {

/// contracted Gaussian geminal = \f$ \sum_i c_i \exp(- \alpha r_{12}^2) \f$,
/// represented as a vector of
/// {\f$ \alpha_i \f$, \f$ c_i \f$ } pairs
typedef std::vector<std::pair<double, double>> ContractedGaussianGeminal;

constexpr size_t num_geometrical_derivatives(size_t ncenter,
                                             size_t deriv_order) {
  return (deriv_order > 0)
             ? (num_geometrical_derivatives(ncenter, deriv_order - 1) *
                (3 * ncenter + deriv_order - 1)) /
                   deriv_order
             : 1;
}

template <typename T, unsigned N>
__libint2_engine_inline
    typename std::remove_all_extents<T>::type* to_ptr1(T (&a)[N]);

/// \brief types of operators (operator sets) supported by Engine.

/// \warning These must start with 0 and appear in same order as elements of
/// BOOST_PP_NBODY_OPERATOR_LIST preprocessor macro (aliases do not need to be
/// included). \warning for the sake of nbody() order operators by # of
/// particles
enum class Operator {
  /// overlap
  overlap = 0,
  /// electronic kinetic energy, i.e. \f$ -\frac{1}{2} \nabla^2 \f$
  kinetic,
  /// Coulomb potential due to point charges
  nuclear,
  /// erf-attenuated point-charge Coulomb operator,
  /// \f$ \mathrm{erf}(\omega r)/r \f$
  erf_nuclear,
  /// erfc-attenuated point-charge Coulomb operator,
  /// \f$ \mathrm{erfc}(\omega r)/r \f$
  erfc_nuclear,
  /// mixture of erf_nuclear and erfc_nuclear
  /// \f$ (\lambda \mathrm{erf}(\omega r) + \sigma \mathrm{erfc}(\omega r) )/r
  /// \f$
  erfx_nuclear,
  //! overlap + (Cartesian) electric dipole moment,
  //! \f$ x_O, y_O, z_O \f$, where
  //! \f$ x_O \equiv x - O_x \f$ is relative to
  //! origin \f$ \vec{O} \f$
  emultipole1,
  //! emultipole1 + (Cartesian) electric quadrupole moment,
  //! \f$ x^2, xy, xz, y^2, yz, z^2 \f$
  emultipole2,
  //! emultipole2 + (Cartesian) electric octupole moment,
  //! \f$ x^3, x^2y, x^2z, xy^2, xyz, xz^2, y^3, y^2z, yz^2, z^3 \f$
  emultipole3,
  //! (electric) spherical multipole moments,
  //! \f$ O_{l,m} \equiv \mathcal{N}^{\text{sign}(m)}_{l,|m|} \f$ where \f$ \f$
  //! \mathcal{N}^{\pm}_{l,m} \f$
  //! is defined in J.M. Pérez-Jordá and W. Yang, J Chem Phys 104, 8003 (1996),
  //! DOI 10.1063/1.468354 .
  //! To obtain the real solid harmonics \f$ C^m_l \f$ and \f$ S^m_l \f$ defined
  //! in https://en.wikipedia.org/wiki/Solid_harmonics
  //! multiply these harmonics by \f$ (-1)^m \sqrt{(2 - \delta_{m,0}) (l + |m|)!
  //! (l - |m|)!} \f$ .
  //! The operator set includes multipoles of order up to \f$ l_{\rm max} = \f$
  //! MULTIPOLE_MAX_ORDER (for a total of \f$ (l_{\rm max}+1)^2 \f$ operators),
  //! in the order of increasing \c l , with the operators of same \c l but
  //! different \c m ordered according to the solid harmonics ordering
  //! CCA standard (see macro FOR_SOLIDHARM_STANDARD in shgshell_ordering.h.in).
  //! For example, the operators will appear in the following order
  //! \f$ \mathcal{N}^+_{0,0} , \mathcal{N}^-_{1,1}, \mathcal{N}^+_{1,0},
  //! \mathcal{N}^+_{1,1}, \mathcal{N}^-_{2,2}, \mathcal{N}^-_{2,1},
  //! \mathcal{N}^+_{2,0}, \mathcal{N}^+_{2,1}, \mathcal{N}^+_{2,2}. \dots \f$ .
  //! Previous to cdbb9f3 released in v2.8.0, Standard -or- Gaussian ordering
  //! could be be specified at generator/compiler configure time.
  sphemultipole,
  /// The four components of σp . V . σp, where V is the nuclear potential.
  opVop,
  /// \f$ \delta(\vec{r}_1 - \vec{r}_2) \f$
  delta,
  /// (2-body) Coulomb operator = \f$ r_{12}^{-1} \f$
  coulomb,
  /// alias for Operator::coulomb
  r12_m1 = coulomb,
  /// contracted Gaussian geminal
  cgtg,
  /// contracted Gaussian geminal times Coulomb
  cgtg_x_coulomb,
  /// |Delta . contracted Gaussian geminal|^2
  delcgtg2,
  /// anti-Coulomb operator, \f$ r_{12} \f$
  r12,
  /// alias for Operator::r12
  r12_1 = r12,
  /// erf-attenuated Coulomb operator,
  /// \f$ \mathrm{erf}(\omega r)/r \f$
  erf_coulomb,
  /// erfc-attenuated Coulomb operator,
  /// \f$ \mathrm{erfc}(\omega r)/r \f$
  erfc_coulomb,
  /// mixture of erf_coulomb and erfc_coulomb
  /// \f$ (\lambda \mathrm{erf}(\omega r) + \sigma \mathrm{erfc}(\omega r) )/r
  /// \f$
  erfx_coulomb,
  /// Slater-type geminal, \f$ \mathrm{exp}(-\zeta r_{12}) \f$
  stg,
  /// Slater-type geminal times Coulomb, , \f$ \mathrm{exp}(-\zeta r_{12}) /
  /// r_{12} \f$
  stg_x_coulomb,
  /// alias for Operator::stg_x_coulomb
  yukawa = stg_x_coulomb,
  // do not modify this
  invalid = -1,
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!keep this
  // updated!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  first_1body_oper = overlap,
  last_1body_oper = opVop,
  first_2body_oper = delta,
  last_2body_oper = stg_x_coulomb,
  first_oper = first_1body_oper,
  last_oper = last_2body_oper
};

/// @param[in] oper an Operator object
/// @return the particle rank of \c oper
/// @throw std::logic_error if an invalid @c oper given
inline constexpr int rank(Operator oper) {
  return (oper >= Operator::first_1body_oper &&
          oper <= Operator::last_1body_oper)
             ? 1
             : ((oper >= Operator::first_2body_oper &&
                 oper <= Operator::last_2body_oper)
                    ? 2
                    : throw std::logic_error(
                          "rank(Operator): invalid operator given"));
}

namespace detail {
struct default_operator_traits {
  typedef struct {
  } oper_params_type;
  static oper_params_type default_params() { return oper_params_type{}; }
  static constexpr auto nopers = 1u;
  // N.B.: Below field means we *should* template specialize operator_traits for
  // Operator::kinetic, but L2 doesn't use that anywhere.
  static constexpr auto intrinsic_deriv_order = 0u;
  struct _core_eval_type {
    template <typename... params>
    static std::shared_ptr<const _core_eval_type> instance(params...) {
      return nullptr;
    }
  };
  using core_eval_type = const _core_eval_type;
};
}  // namespace detail

/// describes operator set \c Op
/// @tparam Op a value of type Operator
/// \note default describes operator set of size 1 that takes trivial \c
/// oper_params_type and \c core_eval_type;
///       needs to be specialized for some operator types
template <Operator Op>
struct operator_traits : public detail::default_operator_traits {};

template <>
struct operator_traits<Operator::nuclear>
    : public detail::default_operator_traits {
  /// point charges and their positions
  typedef std::vector<std::pair<scalar_type, std::array<scalar_type, 3>>>
      oper_params_type;
  static oper_params_type default_params() { return oper_params_type{}; }
#ifndef LIBINT_USER_DEFINED_REAL
  typedef const libint2::FmEval_Chebyshev7<scalar_type> core_eval_type;
#else
  typedef const libint2::FmEval_Reference<scalar_type> core_eval_type;
#endif
};
template <>
struct operator_traits<Operator::opVop>
    : public operator_traits<Operator::nuclear> {
  static constexpr auto nopers = 4;
  static constexpr auto intrinsic_deriv_order = 2;
};

template <>
struct operator_traits<Operator::erf_nuclear>
    : public detail::default_operator_traits {
  /// the attenuation parameter (0 = zero potential, +infinity = no attenuation)
  /// + point charges and positions
  typedef std::tuple<scalar_type, typename operator_traits<
                                      Operator::nuclear>::oper_params_type>
      oper_params_type;
  static oper_params_type default_params() {
    return std::make_tuple(
        0, operator_traits<Operator::nuclear>::default_params());
  }
  typedef const libint2::GenericGmEval<
      libint2::os_core_ints::erf_coulomb_gm_eval<scalar_type>>
      core_eval_type;
};

template <>
struct operator_traits<Operator::erfc_nuclear>
    : public detail::default_operator_traits {
  /// the attenuation parameter (0 = no attenuation, +infinity = zero potential)
  /// + point charges and positions
  typedef typename operator_traits<Operator::erf_nuclear>::oper_params_type
      oper_params_type;
  static oper_params_type default_params() {
    return std::make_tuple(
        0, operator_traits<Operator::nuclear>::default_params());
  }
  typedef const libint2::GenericGmEval<
      libint2::os_core_ints::erfx_coulomb_gm_eval<scalar_type>>
      core_eval_type;
};

template <>
struct operator_traits<Operator::erfx_nuclear>
    : public detail::default_operator_traits {
  /// {{attenuation parameter (inverse lengthscale) + coefficient of erf
  /// (long-range) + coefficient of erfc (short-range)}, nuclear
  /// charges/positions}
  typedef std::tuple<
      std::array<scalar_type, 3>,
      typename operator_traits<Operator::nuclear>::oper_params_type>
      oper_params_type;
  static oper_params_type default_params() {
    return std::make_tuple(
        std::array<scalar_type, 3>{0., 0., 0.},
        operator_traits<Operator::nuclear>::default_params());
  }
  typedef const libint2::GenericGmEval<
      libint2::os_core_ints::erfx_coulomb_gm_eval<scalar_type>>
      core_eval_type;
};

template <>
struct operator_traits<Operator::emultipole1>
    : public detail::default_operator_traits {
  /// Cartesian coordinates of the origin with respect to which the dipole
  /// moment is defined
  typedef std::array<double, 3> oper_params_type;
  static oper_params_type default_params() {
    return oper_params_type{{0, 0, 0}};
  }
  static constexpr auto nopers = 4u;  //!< overlap + 3 dipole components
};
template <>
struct operator_traits<Operator::emultipole2>
    : public operator_traits<Operator::emultipole1> {
  static constexpr auto nopers =
      operator_traits<Operator::emultipole1>::nopers +
      6;  //!< overlap + 3 dipoles + 6 quadrupoles
};
template <>
struct operator_traits<Operator::emultipole3>
    : public operator_traits<Operator::emultipole1> {
  static constexpr auto nopers =
      operator_traits<Operator::emultipole2>::nopers + 10;
};
template <>
struct operator_traits<Operator::sphemultipole>
    : public operator_traits<Operator::emultipole1> {
  static constexpr auto nopers =
      (MULTIPOLE_MAX_ORDER + 1) * (MULTIPOLE_MAX_ORDER + 1);
};

template <>
struct operator_traits<Operator::coulomb>
    : public detail::default_operator_traits {
#ifndef LIBINT_USER_DEFINED_REAL
  typedef const libint2::FmEval_Chebyshev7<scalar_type> core_eval_type;
#else
  typedef const libint2::FmEval_Reference<scalar_type> core_eval_type;
#endif
};
namespace detail {
template <int K>
struct cgtg_operator_traits : public detail::default_operator_traits {
  typedef libint2::GaussianGmEval<scalar_type, K> core_eval_type;
  typedef ContractedGaussianGeminal oper_params_type;
};
}  // namespace detail
template <>
struct operator_traits<Operator::cgtg>
    : public detail::cgtg_operator_traits<0> {};
template <>
struct operator_traits<Operator::cgtg_x_coulomb>
    : public detail::cgtg_operator_traits<-1> {};
template <>
struct operator_traits<Operator::delcgtg2>
    : public detail::cgtg_operator_traits<2> {};

template <>
struct operator_traits<Operator::delta>
    : public detail::default_operator_traits {
  typedef const libint2::GenericGmEval<
      libint2::os_core_ints::delta_gm_eval<scalar_type>>
      core_eval_type;
};

template <>
struct operator_traits<Operator::r12> : public detail::default_operator_traits {
  typedef const libint2::GenericGmEval<
      libint2::os_core_ints::r12_xx_K_gm_eval<scalar_type, 1>>
      core_eval_type;
};

template <>
struct operator_traits<Operator::erf_coulomb>
    : public detail::default_operator_traits {
  /// the attenuation parameter (inverse lengthscale, 0 = zero potential,
  /// +infinity = no attenuation)
  typedef scalar_type oper_params_type;
  static oper_params_type default_params() { return oper_params_type{0}; }
  typedef const libint2::GenericGmEval<
      libint2::os_core_ints::erf_coulomb_gm_eval<scalar_type>>
      core_eval_type;
};
template <>
struct operator_traits<Operator::erfc_coulomb>
    : public detail::default_operator_traits {
  /// the attenuation parameter (inverse lengthscale, 0 = no attenuation,
  /// +infinity = zero potential)
  typedef scalar_type oper_params_type;
  static oper_params_type default_params() { return oper_params_type{0}; }
  typedef const libint2::GenericGmEval<
      libint2::os_core_ints::erfx_coulomb_gm_eval<scalar_type>>
      core_eval_type;
};
template <>
struct operator_traits<Operator::erfx_coulomb>
    : public detail::default_operator_traits {
  /// the attenuation parameter (inverse lengthscale) + coefficient of erf
  /// (long-range) + coefficient of erfc (short-range)
  typedef std::array<scalar_type, 3> oper_params_type;
  static oper_params_type default_params() {
    return oper_params_type{0., 0., 0.};
  }
  typedef const libint2::GenericGmEval<
      libint2::os_core_ints::erfx_coulomb_gm_eval<scalar_type>>
      core_eval_type;
};

template <>
struct operator_traits<Operator::stg> : public detail::default_operator_traits {
  /// the attenuation parameter (0 = constant, infinity = delta-function)
  typedef scalar_type oper_params_type;
  static oper_params_type default_params() { return oper_params_type{0}; }
  typedef const libint2::TennoGmEval<scalar_type> core_eval_type;
};

template <>
struct operator_traits<Operator::stg_x_coulomb>
    : public detail::default_operator_traits {
  /// the attenuation parameter (0 = coulomb, infinity = delta-function)
  typedef scalar_type oper_params_type;
  static oper_params_type default_params() { return oper_params_type{0}; }
  typedef const libint2::TennoGmEval<scalar_type> core_eval_type;
};

/// the runtime version of \c operator_traits<oper>::default_params()
libint2::any default_params(const Operator& oper);

namespace detail {
template <typename core_eval_type>
using __core_eval_pack_type =
    compressed_pair<std::shared_ptr<core_eval_type>,
                    libint2::detail::CoreEvalScratch<core_eval_type>>;
template <Operator Op>
using core_eval_pack_type =
    __core_eval_pack_type<typename operator_traits<Op>::core_eval_type>;
}  // namespace detail

#define BOOST_PP_NBODY_BRAKET_MAX_INDEX 4

/// @param[in] braket a BraKet object
/// @return rank of @c braket
/// @throw std::logic_error if invalid @c braket given
inline constexpr int rank(BraKet braket) {
  return (braket == BraKet::x_x || braket == BraKet::xs_xs)
             ? 2
             : ((braket == BraKet::xs_xx || braket == BraKet::xx_xs)
                    ? 3
                    : ((braket == BraKet::xx_xx)
                           ? 4
                           : throw std::logic_error(
                                 "rank(BraKet): invalid braket given")));
}

/// @param[in] oper an Operator object
/// @return the default braket for @c oper
/// @throw std::logic_error if invalid @c oper given
inline constexpr BraKet default_braket(Operator oper) {
  return (rank(oper) == 1)
             ? BraKet::x_x
             : ((rank(oper) == 2)
                    ? BraKet::xx_xx
                    : throw std::logic_error(
                          "default_braket(Operator): invalid operator given"));
}

constexpr size_t nopers_2body = static_cast<int>(Operator::last_2body_oper) -
                                static_cast<int>(Operator::first_2body_oper) +
                                1;
constexpr size_t nbrakets_2body = static_cast<int>(BraKet::last_2body_braket) -
                                  static_cast<int>(BraKet::first_2body_braket) +
                                  1;
constexpr size_t nderivorders_2body = LIBINT2_MAX_DERIV_ORDER + 1;

/**
 * Engine computes integrals of operators (or operator sets) specified by
 * combination of Operator and BraKet.
 * This class deprecates OneBodyEngine and TwoBodyEngine.
 */
class Engine {
 private:
  typedef struct {
  } empty_pod;

 public:
  static constexpr auto max_ntargets =
      std::extent<decltype(std::declval<Libint_t>().targets), 0>::value;
  using target_ptr_vec =
      std::vector<const value_type*,
                  detail::ext_stack_allocator<const value_type*, max_ntargets>>;

  /// creates a default Engine that cannot be used for computing integrals;
  /// to be used as placeholder for copying a usable engine, OR for cleanup of
  /// thread-local data. Nontrivial use of a default-initialized Engine will
  /// cause an exception of type Engine::using_default_initialized
  Engine()
      : oper_(Operator::invalid),
        braket_(BraKet::invalid),
        primdata_(),
        stack_size_(0),
        lmax_(-1),
        deriv_order_(0),
        cartesian_shell_normalization_(CartesianShellNormalization::standard),
        scale_(1) {
    set_precision(std::numeric_limits<scalar_type>::epsilon());
  }

  // clang-format off
  /// Constructs a (usable) Engine

  /// \param oper a value of Operator type
  /// \param max_nprim the maximum number of primitives per contracted Gaussian
  /// shell (must be greater than 0)
  /// \param max_l the maximum angular momentum of Gaussian shell (>=0)
  /// \throw Engine::lmax_exceeded if \c max_l exceeds the angular momentum
  /// limit of the library
  /// \param deriv_order if not 0, will compute geometric derivatives of
  ///        Gaussian integrals of order \c deriv_order
  ///        (default=0, i.e. compute nondifferentiated integrals)
  /// \param precision specifies the target precision with which the integrals
  ///        will be computed; the default is the
  ///        `std::numeric_limits<scalar_type>::epsilon()`.
  ///        Currently the precision control is implemented
  ///        for two-body integrals only.
  ///        \sa Engine::set_precision()
  /// \param params a value of type `Engine::operator_traits<oper>::oper_params_type`
  ///               specifying the parameters of
  ///               the operator set, e.g. position and magnitude of the charges
  ///               creating the Coulomb potential
  ///               for `oper == Operator::nuclear`, etc.
  ///               For most values of \c oper
  ///               this is not needed.
  ///               \sa Engine::operator_traits
  /// \param braket a value of BraKet type
  /// \warning currently only one-contraction Shell objects are supported; i.e.
  /// generally-contracted Shells are not yet supported
  // clang-format on
  template <typename Params = empty_pod>
  Engine(Operator oper, size_t max_nprim, int max_l, int deriv_order = 0,
         scalar_type precision = std::numeric_limits<scalar_type>::epsilon(),
         Params params = empty_pod(), BraKet braket = BraKet::invalid,
         ScreeningMethod screening_method = default_screening_method())
      : oper_(oper),
        braket_(braket),
        primdata_(),
        spbra_(max_nprim),
        spket_(max_nprim),
        stack_size_(0),
        lmax_(max_l),
        deriv_order_(deriv_order),
        screening_method_(screening_method),
        cartesian_shell_normalization_(CartesianShellNormalization::standard),
        scale_(1),
        params_(enforce_params_type(oper, params)) {
    set_precision(precision);
    assert(max_nprim > 0);
    initialize(max_nprim);
    core_eval_pack_ = make_core_eval_pack(oper);  // must follow initialize() to
                                                  // ensure default braket_ has
                                                  // been set
    init_core_ints_params(params_);
  }

  /// move constructor
  // intel does not accept "move ctor = default"
  Engine(Engine&& other)
      : oper_(other.oper_),
        braket_(other.braket_),
        primdata_(std::move(other.primdata_)),
        spbra_(std::move(other.spbra_)),
        spket_(std::move(other.spket_)),
        stack_size_(other.stack_size_),
        lmax_(other.lmax_),
        hard_lmax_(other.hard_lmax_),
        hard_default_lmax_(other.hard_default_lmax_),
        deriv_order_(other.deriv_order_),
        precision_(other.precision_),
        ln_precision_(other.ln_precision_),
        screening_method_(other.screening_method_),
        cartesian_shell_normalization_(other.cartesian_shell_normalization_),
        scale_(other.scale_),
        core_eval_pack_(std::move(other.core_eval_pack_)),
        params_(std::move(other.params_)),
        core_ints_params_(std::move(other.core_ints_params_)),
        targets_(std::move(other.targets_)),
        set_targets_(other.set_targets_),
        scratch_(std::move(other.scratch_)),
        scratch2_(other.scratch2_),
        buildfnptrs_(other.buildfnptrs_) {
    // leave other in an unusable (but valid) state
    other.oper_ = Operator::invalid;
    other.braket_ = BraKet::invalid;
    other.scratch2_ = nullptr;
  }

  /// (deep) copy constructor
  Engine(const Engine& other)
      : oper_(other.oper_),
        braket_(other.braket_),
        primdata_(other.primdata_.size()),
        spbra_(other.spbra_),
        spket_(other.spket_),
        stack_size_(other.stack_size_),
        lmax_(other.lmax_),
        deriv_order_(other.deriv_order_),
        precision_(other.precision_),
        ln_precision_(other.ln_precision_),
        screening_method_(other.screening_method_),
        cartesian_shell_normalization_(other.cartesian_shell_normalization_),
        scale_(other.scale_),
        core_eval_pack_(other.core_eval_pack_),
        params_(other.params_),
        core_ints_params_(other.core_ints_params_) {
    initialize();
  }

  ~Engine() { finalize(); }

  /// move assignment is default
  Engine& operator=(Engine&& other) {
    oper_ = other.oper_;
    braket_ = other.braket_;
    primdata_ = std::move(other.primdata_);
    spbra_ = std::move(other.spbra_);
    spket_ = std::move(other.spket_);
    stack_size_ = other.stack_size_;
    lmax_ = other.lmax_;
    hard_lmax_ = other.hard_lmax_;
    hard_default_lmax_ = other.hard_default_lmax_;
    deriv_order_ = other.deriv_order_;
    precision_ = other.precision_;
    ln_precision_ = other.ln_precision_;
    screening_method_ = other.screening_method_;
    cartesian_shell_normalization_ = other.cartesian_shell_normalization_;
    scale_ = other.scale_;
    core_eval_pack_ = std::move(other.core_eval_pack_);
    params_ = std::move(other.params_);
    core_ints_params_ = std::move(other.core_ints_params_);
    targets_ = std::move(other.targets_);
    set_targets_ = other.set_targets_;
    scratch_ = std::move(other.scratch_);
    scratch2_ = other.scratch2_;
    buildfnptrs_ = other.buildfnptrs_;
    // leave other in an unusable state
    other.oper_ = Operator::invalid;
    other.braket_ = BraKet::invalid;
    other.scratch2_ = nullptr;
    return *this;
  }

  /// (deep) copy assignment
  Engine& operator=(const Engine& other) {
    oper_ = other.oper_;
    braket_ = other.braket_;
    primdata_.resize(other.primdata_.size());
    spbra_ = other.spbra_;
    spket_ = other.spket_;
    stack_size_ = other.stack_size_;
    lmax_ = other.lmax_;
    deriv_order_ = other.deriv_order_;
    precision_ = other.precision_;
    ln_precision_ = other.ln_precision_;
    screening_method_ = other.screening_method_;
    cartesian_shell_normalization_ = other.cartesian_shell_normalization_;
    scale_ = other.scale_;
    core_eval_pack_ = other.core_eval_pack_;
    params_ = other.params_;
    core_ints_params_ = other.core_ints_params_;
    initialize();
    return *this;
  }

  /// @return the particle rank of the operator
  int operator_rank() const { return rank(oper_); }

  /// @return rank of the braket (e.g., 2 for (a|b), 3 for (a|bc), etc.)
  int braket_rank() const { return rank(braket_); }

  /// @return the operator
  Operator oper() const { return oper_; }

  /// @return the braket
  BraKet braket() const { return braket_; }

  /// @return the order of geometrical derivatives
  int deriv_order() const { return deriv_order_; }

  /// (re)sets operator type to @c new_oper
  /// @param[in] new_oper Operator whose integrals will be computed with the
  /// next call to Engine::compute()
  /// @note this resets braket and params to their respective defaults for @c
  /// new_oper
  Engine& set(Operator new_oper) {
    if (oper_ != new_oper) {
      oper_ = new_oper;
      braket_ = default_braket(oper_);
      params_ = default_params(oper_);
      initialize();
      core_eval_pack_ = make_core_eval_pack(oper_);  // must follow initialize()
                                                     // to ensure default
                                                     // braket_ has been set
    }
    return *this;
  }

  /// (re)sets braket type to @c new_braket
  /// @param[in] new_braket integrals BraKet that will be computed with the next
  /// call to Engine::compute()
  Engine& set(BraKet new_braket) {
    if (braket_ != new_braket) {
      braket_ = new_braket;
      initialize();
    }
    return *this;
  }

  /// resets operator parameters; this may be useful e.g. if need to compute
  /// Coulomb potential
  /// integrals over batches of charges for the sake of parallelism.
  template <typename Params>
  Engine& set_params(const Params& params) {
    params_ = params;
    init_core_ints_params(params_);
    reset_scratch();
    return *this;
  }

  /// @return the maximum number of primitives that this engine can handle
  std::size_t max_nprim() const {
    assert(spbra_.primpairs.size() == spket_.primpairs.size());
    return static_cast<std::size_t>(std::sqrt(spbra_.primpairs.size()));
  }

  /// reset the maximum number of primitives
  /// @param[in] n the maximum number of primitives
  /// @note left unchanged if the value returned by Engine::max_nprim is greater
  /// than @c n
  Engine& set_max_nprim(std::size_t n) {
    if (n * n > spbra_.primpairs.size()) {
      spbra_.resize(n);
      spket_.resize(n);
      initialize(n);
    }
    return *this;
  }

  /// @return the maximum angular momentum that this engine can handle
  std::size_t max_l() const { return lmax_; }

  /// reset the maximum angular momentum
  /// @param[in] L the maximum angular momentum
  /// @note left unchanged if the value returned by Engine::max_l is greater
  /// than @c L
  Engine& set_max_l(std::size_t L) {
    if (L >= static_cast<std::size_t>(lmax_)) {
      lmax_ = L;
      initialize();
    }
    return *this;
  }

  /// returns a vector that will hold pointers to shell sets computed with
  /// Engine::compute()
  /// or other compute functions. Only need to get this vector once, but the
  /// values will change
  /// after every compute() call.
  const target_ptr_vec& results() const { return targets_; }

  /// reports the number of shell sets that each call to compute() produces.
  /// this depends on the order of geometrical derivatives requested and
  /// on the operator set. \sa compute_nshellsets()
  unsigned int nshellsets() const { return targets_.size(); }

  /// Computes target shell sets of integrals.

  /// @return vector of pointers to target shell sets, the number of sets =
  /// Engine::nshellsets();
  ///         if the first pointer equals \c nullptr then all elements were
  ///         screened out.
  /// \note resulting shell sets are stored in row-major order.
  /// \note Call Engine::compute1() or Engine::compute2() directly to avoid
  /// extra copies.
  template <typename... ShellPack>
  __libint2_engine_inline const target_ptr_vec& compute(
      const libint2::Shell& first_shell, const ShellPack&... rest_of_shells);

  /// Computes target shell sets of 1-body integrals.
  /// @param[in] s1 bra shell
  /// @param[in] s2 ket shell
  /// @return vector of pointers to target shell sets, the number of sets =
  /// Engine::nshellsets() \note resulting shell sets are stored in row-major
  /// order
  __libint2_engine_inline const target_ptr_vec& compute1(
      const libint2::Shell& s1, const libint2::Shell& s2);

  /// Computes target shell sets of 2-body integrals, @code  @endcode
  /// @note result is stored in the "chemists"/Mulliken form, @code (bra1
  /// bra2|ket1 ket2) @endcode; the "physicists"/Dirac bra-ket form is @code
  /// <bra1 ket1|bra2 ket2> @endcode, where @c bra1 and @c bra2 refer to
  /// particle 1, and @c ket1 and @c ket2 refer to particle 2.
  /// @tparam oper operator
  /// @tparam braket the integral type
  /// @tparam deriv_order the derivative order, values greater than 2 not yet
  /// supported
  /// @param[in] bra1 the first shell in the Mulliken bra
  /// @param[in] bra2 the second shell in the Mulliken bra
  /// @param[in] ket1 the first shell in the Mulliken ket
  /// @param[in] ket2 the second shell in the Mulliken ket
  /// @param[in] spbra ShellPair data for shell pair @c {bra1,bra2}
  /// @param[in] spket ShellPair data for shell pair @c {ket1,ket2}
  /// @return vector of pointers to target shell sets, the number of sets =
  /// Engine::nshellsets();
  ///         if the first pointer equals \c nullptr then all elements were
  ///         screened out.
  /// @note the result integrals are packed in shell sets in row-major order
  /// (i.e. function index of shell ket2 has stride 1, function index of shell
  /// ket1 has stride nket2, etc.).
  /// @note internally the integrals are evaluated with shells permuted to
  /// according to the canonical order predefined at the library generation time
  /// (see macro @c LIBINT_SHELL_SET ). To minimize the overhead it is
  /// recommended to organize your shell loop nests in the order best suited for
  /// your particular instance of Libint library.
  template <Operator oper, BraKet braket, size_t deriv_order>
  __libint2_engine_inline const target_ptr_vec& compute2(
      const Shell& bra1, const Shell& bra2, const Shell& ket1,
      const Shell& ket2, const ShellPair* spbra = nullptr,
      const ShellPair* spket = nullptr);

  typedef const target_ptr_vec& (Engine::*compute2_ptr_type)(
      const Shell& bra1, const Shell& bra2, const Shell& ket1,
      const Shell& ket2, const ShellPair* spbra, const ShellPair* spket);

  // clang-format off
  /** this specifies target precision for computing the integrals, i.e.
   *  the target absolute (i.e., not relative) error of the integrals.
   *  It is used to screen out primitive integrals. For some screening
   *  methods precision can be almost guaranteed (due to finite precision
   *  of the precomputed interpolation tables used to evaluate the core integrals
   *  it is not in general possible to guarantee precision rigorously).
   *
   *  @param[in] prec the target precision
   *  @sa ScreeningMethod
   */
  // clang-format on
  Engine& set_precision(scalar_type prec) {
    if (prec <= 0.) {
      precision_ = 0.;
      ln_precision_ = std::numeric_limits<scalar_type>::lowest();
      return *this;
    }
    if (prec > 0.) {  // split single if/else to avoid speculative execution of
                      // else branch on Apple M1 with apple clang 13.0.0
      precision_ = prec;
      using std::log;
      ln_precision_ = log(precision_);
      return *this;
    }
    abort();
  }
  /// @return the target precision for computing the integrals
  /// @sa set_precision(scalar_type)
  scalar_type precision() const { return precision_; }

  /// @param screening_method method used to screen primitive contributions
  Engine& set(ScreeningMethod screening_method) {
    screening_method_ = screening_method;
    return *this;
  }

  /// @return the method used to screen primitive contributions
  ScreeningMethod screening_method() const { return screening_method_; }

  /// @return the Cartesian Gaussian normalization convention
  CartesianShellNormalization cartesian_shell_normalization() const {
    return cartesian_shell_normalization_;
  }

  /// @param[in] norm the normalization convention for Cartesian Gaussians
  /// @note the default convention is CartesianShellNormalization_Standard
  /// @return reference to @c this for daisy-chaining
  Engine& set(CartesianShellNormalization norm) {
    cartesian_shell_normalization_ = norm;
    return *this;
  }

  /// @return the scaling factor by which the target integrals are multiplied
  scalar_type prescaled_by() const { return scale_; }

  /// @param[in] sc the scaling factor by which the target integrals are
  /// multiplied
  /// @note the default factor is 1
  /// @return reference to @c this for daisy-chaining
  Engine& prescale_by(scalar_type sc) {
    scale_ = sc;
    return *this;
  }

  /// prints the contents of timers to @c os
  void print_timers(std::ostream& os = std::cout) {
#ifdef LIBINT2_ENGINE_TIMERS
    os << "timers: prereq = " << timers.read(0);
#ifndef LIBINT2_PROFILE  // if libint's profiling was on, engine's build timer
                         // will include its overhead
    // do not report it, detailed profiling from libint will be printed below
    os << " build = " << timers.read(1);
#endif
    os << " tform = " << timers.read(2) << std::endl;
#endif
#ifdef LIBINT2_PROFILE
    os << "build timers: hrr = " << primdata_[0].timers->read(0)
       << " vrr = " << primdata_[0].timers->read(1) << std::endl;
#endif
#ifdef LIBINT2_ENGINE_TIMERS
#ifdef LIBINT2_ENGINE_PROFILE_CLASS
    for (const auto& p : class_profiles) {
      char buf[1024];
      std::snprintf(buf, sizeof(buf),
                    "{\"%s\", %10.5lf, %10.5lf, %10.5lf, %10.5lf, %ld, %ld},",
                    p.first.to_string().c_str(), p.second.prereqs,
                    p.second.build_vrr, p.second.build_hrr, p.second.tform,
                    p.second.nshellset, p.second.nprimset);
      os << buf << std::endl;
    }
#endif
#endif
  }

  /// Exception class to be used when the angular momentum limit is exceeded.
  class lmax_exceeded : virtual public std::logic_error {
   public:
    lmax_exceeded(const char* task_name, size_t lmax_limit,
                  size_t lmax_requested)
        : std::logic_error(
              "Engine::lmax_exceeded -- angular momentum limit exceeded"),
          lmax_limit_(lmax_limit),
          lmax_requested_(lmax_requested) {
      strncpy(task_name_, task_name, 64);
      task_name_[64] = '\0';
    }
    ~lmax_exceeded() noexcept {}

    const char* task_name() const {
      return static_cast<const char*>(task_name_);
    }
    size_t lmax_limit() const { return lmax_limit_; }
    size_t lmax_requested() const { return lmax_requested_; }

   private:
    char task_name_[65];
    size_t lmax_limit_;
    size_t lmax_requested_;
  };

  /// Exception class to be used when using default-initialized Engine
  class using_default_initialized : virtual public std::logic_error {
   public:
    using_default_initialized()
        : std::logic_error(
              "Engine::using_default_initialized -- attempt to use a "
              "default-initialized Engine") {}
    ~using_default_initialized() noexcept {}
  };

 private:
  Operator oper_;
  BraKet braket_;
  std::vector<Libint_t> primdata_;
  ShellPair spbra_, spket_;
  size_t stack_size_;  // amount allocated by libint2_init_xxx in
                       // primdata_[0].stack
  int lmax_;
  int hard_lmax_;  // max L supported by library for this operator type (i.e.
                   // LIBINT2_MAX_AM_<task><deriv_order>) + 1
  int hard_default_lmax_;  // max L supported by library for default operator
                           // type (i.e. LIBINT2_MAX_AM_default<deriv_order>) +
                           // 1 this is only set and used if
                           // LIBINT2_CENTER_DEPENDENT_MAX_AM_<task><deriv_order>
                           // == 1
  int deriv_order_;
  scalar_type precision_;
  scalar_type ln_precision_;
  ScreeningMethod screening_method_ = ScreeningMethod::Invalid;

  // specifies the normalization convention for Cartesian Gaussians
  CartesianShellNormalization cartesian_shell_normalization_;

  // the target integrals are scaled by this factor (of course it is actually
  // the core integrals that are scaled)
  scalar_type scale_;

  any core_eval_pack_;

  /// operator params
  any params_;  // operator params
  /// for some operators need core ints params that are computed from operator
  /// params,
  /// e.g. integrals of \f$ f_{12}^2 \f$ are computed from parameters of \f$
  /// f_{12} \f$
  any core_ints_params_;
  /// makes core ints params from the operator params
  void init_core_ints_params(const any& params);

  /// pointers to target shell sets, size is updated by reset_scratch()
  /// targets_.size() is returned by nshellsets()
  target_ptr_vec targets_;
  /// true if targets_ does not point primdata_[0].targets
  /// hence must set its contents explicitly
  bool set_targets_;

  std::vector<value_type>
      scratch_;  // for transposes and/or transforming to solid harmonics
  value_type* scratch2_;  // &scratch_[0] points to the first block large enough
                          // to hold all target ints
  // scratch2_ points to second such block. It could point into scratch_ or at
  // primdata_[0].stack
  typedef void (*buildfnptr_t)(const Libint_t*);
  buildfnptr_t* buildfnptrs_;

  /// reports the number of shell sets that each call to compute() produces.
  unsigned int compute_nshellsets() const {
    const unsigned int num_operator_geometrical_derivatives =
        (oper_ == Operator::nuclear || oper_ == Operator::erf_nuclear ||
         oper_ == Operator::erfc_nuclear || oper_ == Operator::erfx_nuclear)
            ? this->nparams()
            : 0;
    const auto ncenters = braket_rank() + num_operator_geometrical_derivatives;
    return nopers() * num_geometrical_derivatives(ncenters, deriv_order_);
  }

  void reset_scratch() {
    const auto nshsets = compute_nshellsets();
    targets_.resize(nshsets);
    set_targets_ =
        (&targets_[0] != const_cast<const value_type**>(primdata_[0].targets));
    const auto ncart_max = (lmax_ + 1) * (lmax_ + 2) / 2;
    const auto target_shellset_size =
        nshsets * std::pow(ncart_max, braket_rank());
    // need to be able to hold 2 sets of target shellsets: the worst case occurs
    // when dealing with
    // 1-body Coulomb ints derivatives ... have 2+natom 1st-order derivative
    // sets (and many more of 2nd and higher) that are stored in scratch then
    // need to transform to solids. To avoid copying back and forth make sure
    // that there is enough room to transform all ints and save them in correct
    // order in single pass
    const auto need_extra_large_scratch = stack_size_ < target_shellset_size;
    scratch_.resize(need_extra_large_scratch ? 2 * target_shellset_size
                                             : target_shellset_size);
    scratch2_ = need_extra_large_scratch ? &scratch_[target_shellset_size]
                                         : primdata_[0].stack;
  }

  __libint2_engine_inline void compute_primdata(Libint_t& primdata,
                                                const Shell& s1,
                                                const Shell& s2, size_t p1,
                                                size_t p2, size_t oset);

 public:
  /// 3-dim array of pointers to help dispatch efficiently based on oper_,
  /// braket_, and deriv_order_
  /// @note public since 2.7.0 to support efficient dispatch in user code
  __libint2_engine_inline const std::vector<Engine::compute2_ptr_type>&
  compute2_ptrs() const;

 private:
  // max_nprim=0 avoids resizing primdata_
  __libint2_engine_inline void initialize(size_t max_nprim = 0);
  // generic _initializer
  __libint2_engine_inline void _initialize();

  void finalize() {
    if (primdata_.size() != 0) {
      libint2_cleanup_default(&primdata_[0]);
    }
  }  // finalize()

  //-------
  // utils
  //-------
  unsigned int nparams() const;
  unsigned int nopers() const;
  unsigned int intrinsic_deriv_order() const;
  /// if Params == operator_traits<oper>::oper_params_type, will return
  /// \c any(params)
  /// else will set return \c any initialized with default value for
  /// \c operator_traits<type>::oper_params_type
  /// @param throw_if_wrong_type if true, and Params !=
  /// operator_traits<type>::oper_params_type, will throw std::bad_cast
  template <typename Params>
  static any enforce_params_type(
      Operator oper, const Params& params,
      bool throw_if_wrong_type = !std::is_same<Params, empty_pod>::value);
  /// @return core eval pack corresponding to
  /// operator_traits<oper>::core_eval_type
  any make_core_eval_pack(Operator oper) const;

  //-------
  // profiling
  //-------
  static const bool skip_core_ints = false;
};  // struct Engine

}  // namespace libint2

#ifndef LIBINT2_DOES_NOT_INLINE_ENGINE
#include "./engine.impl.h"
#endif

#endif /* _libint2_src_lib_libint_engine_h_ */
