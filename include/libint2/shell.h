/*
 *  Copyright (C) 2004-2023 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_lib_libint_shell_h_
#define _libint2_src_lib_libint_shell_h_

#include <libint2/util/cxxstd.h>
#if LIBINT2_CPLUSPLUS_STD < 2011
#error "libint2/shell.h requires C++11 support"
#endif

#include <libint2.h>
#include <libint2/util/small_vector.h>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

namespace libint2 {

namespace math {
/// fac[k] = k!
static constexpr std::array<int64_t, 21> fac = {{1LL,
                                                 1LL,
                                                 2LL,
                                                 6LL,
                                                 24LL,
                                                 120LL,
                                                 720LL,
                                                 5040LL,
                                                 40320LL,
                                                 362880LL,
                                                 3628800LL,
                                                 39916800LL,
                                                 479001600LL,
                                                 6227020800LL,
                                                 87178291200LL,
                                                 1307674368000LL,
                                                 20922789888000LL,
                                                 355687428096000LL,
                                                 6402373705728000LL,
                                                 121645100408832000LL,
                                                 2432902008176640000LL}};
/// df_Kminus1[k] = (k-1)!!
static constexpr std::array<int64_t, 31> df_Kminus1 = {{1LL,
                                                        1LL,
                                                        1LL,
                                                        2LL,
                                                        3LL,
                                                        8LL,
                                                        15LL,
                                                        48LL,
                                                        105LL,
                                                        384LL,
                                                        945LL,
                                                        3840LL,
                                                        10395LL,
                                                        46080LL,
                                                        135135LL,
                                                        645120LL,
                                                        2027025LL,
                                                        10321920LL,
                                                        34459425LL,
                                                        185794560LL,
                                                        654729075LL,
                                                        3715891200LL,
                                                        13749310575LL,
                                                        81749606400LL,
                                                        316234143225LL,
                                                        1961990553600LL,
                                                        7905853580625LL,
                                                        51011754393600LL,
                                                        213458046676875LL,
                                                        1428329123020800LL,
                                                        6190283353629375LL}};
/// bc(i,j) = binomial coefficient, i! / (j! (i-j)!)
template <typename Int>
int64_t bc(Int i, Int j) {
  assert(i < Int(fac.size()));
  assert(j < Int(fac.size()));
  assert(i >= j);
  return fac[i] / (fac[j] * fac[i - j]);
}
}  // namespace math

/// generally-contracted Solid-Harmonic/Cartesion Gaussian Shell

/** A simple-to-use Gaussian shell. Here's an example of how to create an s+p
 * shell of the STO-3G basis on the oxygen atom located at the origin. \verbatim
 *  auto s = Shell{
 *                  {5.033151300, 1.169596100, 0.380389000},
 *                  {
 *                    {0, false, {-0.09996723, 0.39951283, 0.70011547}},
 *                    {1, false, {0.15591627, 0.60768372, 0.39195739}}
 *                  },
 *                   {{0.0, 0.0, 0.0}}
 *                };
 *  \endverbatim
 *  \note The contraction coefficients correspond to <em>unity-normalized</em>
 * primitives. EMSL Gaussian Basis Set Database, as well as basis set libraries
 * embedded into most quantum chemistry programs use this convention. However,
 * coefficients are automatically converted internally to refer to
 *  normalization-free primitives before computing integrals (see \c
 * Shell::renorm() ).
 */
struct Shell {
  typedef double real_t;

  /// contracted Gaussian = angular momentum + sph/cart flag + contraction
  /// coefficients
  struct Contraction {
    int l;
    bool pure;
    svector<real_t> coeff;

    bool operator==(const Contraction& other) const {
      return &other == this ||
             (l == other.l && pure == other.pure && coeff == other.coeff);
    }
    bool operator!=(const Contraction& other) const {
      return not this->operator==(other);
    }
    size_t cartesian_size() const { return (l + 1) * (l + 2) / 2; }
    size_t size() const { return pure ? (2 * l + 1) : cartesian_size(); }
  };

  svector<real_t> alpha;         //!< exponents
  svector<Contraction> contr;    //!< contractions
  std::array<real_t, 3> O;       //!< origin
  svector<real_t> max_ln_coeff;  //!< maximum ln of (absolute) contraction
                                 //!< coefficient for each primitive

  Shell() = default;
  Shell(const Shell&) = default;
  // intel does not support "move ctor = default"
  Shell(Shell&& other) noexcept
      : alpha(std::move(other.alpha)),
        contr(std::move(other.contr)),
        O(std::move(other.O)),
        max_ln_coeff(std::move(other.max_ln_coeff)) {}
  Shell& operator=(const Shell&) = default;
  // intel does not support "move asgnmt = default"
  Shell& operator=(Shell&& other) noexcept {
    alpha = std::move(other.alpha);
    contr = std::move(other.contr);
    O = std::move(other.O);
    max_ln_coeff = std::move(other.max_ln_coeff);
    return *this;
  }
  /// @param embed_normalization_into_coefficients if true, will embed
  /// normalization factors into coefficients, else will use the coefficients in
  /// @p _contr as given
  Shell(svector<real_t> _alpha, svector<Contraction> _contr,
        std::array<real_t, 3> _O,
        bool embed_normalization_into_coefficients = true)
      : alpha(std::move(_alpha)), contr(std::move(_contr)), O(std::move(_O)) {
    // embed normalization factors into contraction coefficients
    if (embed_normalization_into_coefficients)
      renorm();
    else {
      update_max_ln_coeff();
    }
  }

  Shell& move(std::array<real_t, 3> new_origin) {
    O = std::move(new_origin);
    return *this;
  }

  size_t cartesian_size() const {
    size_t s = 0;
    for (const auto& c : contr) {
      s += c.cartesian_size();
    }
    return s;
  }
  size_t size() const {
    size_t s = 0;
    for (const auto& c : contr) {
      s += c.size();
    }
    return s;
  }

  size_t ncontr() const { return contr.size(); }
  size_t nprim() const { return alpha.size(); }

  bool operator==(const Shell& other) const {
    return &other == this ||
           (O == other.O && alpha == other.alpha && contr == other.contr);
  }
  bool operator!=(const Shell& other) const {
    return not this->operator==(other);
  }

  /// @param l angular momentum quantum number
  /// @return (lower-case) letter symbol corresponding to @p l ; e.g., `s` for
  /// `l=0`, `p` for `l=1`, etc.
  /// @throw std::invalid_argument if \c l is greater than 19
  static char am_symbol(size_t l) {
    static char lsymb[] = "spdfghikmnoqrtuvwxyz";
    assert(l <= 19);
    return lsymb[l];
  }

  /// inverse of am_symbol()
  /// @param am_symbol letter symbol denoting orbital angular momentum @p l ;
  /// e.g., `s` for `l=0`, `p` for `l=1`, etc.
  /// @note this function is case insensitive, i.e. `am_symbol_to_l('s') ==
  /// am_symbol_to_l('S')`
  /// @return angular momentum quantum number
  /// @sa am_symbol()
  static unsigned short am_symbol_to_l(char am_symbol) {
    const char AM_SYMBOL = ::toupper(am_symbol);
    switch (AM_SYMBOL) {
      case 'S':
        return 0;
      case 'P':
        return 1;
      case 'D':
        return 2;
      case 'F':
        return 3;
      case 'G':
        return 4;
      case 'H':
        return 5;
      case 'I':
        return 6;
      case 'K':
        return 7;
      case 'M':
        return 8;
      case 'N':
        return 9;
      case 'O':
        return 10;
      case 'Q':
        return 11;
      case 'R':
        return 12;
      case 'T':
        return 13;
      case 'U':
        return 14;
      case 'V':
        return 15;
      case 'W':
        return 16;
      case 'X':
        return 17;
      case 'Y':
        return 18;
      case 'Z':
        return 19;
      default:
        throw std::invalid_argument{"invalid angular momentum label"};
    }
  }

  struct defaultable_boolean {
    typedef enum { false_value = 0, true_value = 1, default_value = 2 } value_t;
    defaultable_boolean() : value_(default_value) {}
    defaultable_boolean(bool v) : value_(static_cast<value_t>(v ? 1 : 0)) {}
    bool is_default() const { return value_ == default_value; }
    operator bool() const {
      assert(value_ != default_value);
      return value_ == true_value;
    }

   private:
    value_t value_;
  };

  /// sets and/or reports whether the auto-renormalization to unity is set
  /// if called without arguments, returns the current value of the flag
  /// otherwise, will set the flag to \c flag
  /// \note by default, shells WILL be re-normalized to unity
  static bool do_enforce_unit_normalization(
      defaultable_boolean flag = defaultable_boolean()) {
    static bool result{true};
    if (not flag.is_default()) {
      result = flag;
    }
    return result;
  }

  /// @return "unit" Shell, with exponent=0. and coefficient=1., located at the
  /// origin
  static const Shell& unit() {
    static const Shell unitshell{make_unit()};
    return unitshell;
  }

  /// @return the coefficient of primitive \c p in contraction \c c assuming
  /// unit normalized primitive
  ///         (coeff contains coefficients of normalization-free primitives; @sa
  ///         Shell::renorm() )
  real_t coeff_normalized(size_t c, size_t p) const {
    const auto alpha = this->alpha.at(p);
    assert(alpha >= 0.0);
    const auto l = contr.at(c).l;
    assert(l <= 15);  // due to df_Kminus1[] a 64-bit integer type; kinda
                      // ridiculous restriction anyway

    using libint2::math::df_Kminus1;
    const auto sqrt_Pi_cubed = real_t{5.56832799683170784528481798212};
    const auto two_alpha = 2 * alpha;
    const auto two_alpha_to_am32 = pow(two_alpha, l + 1) * sqrt(two_alpha);
    const auto one_over_N = sqrt((sqrt_Pi_cubed * df_Kminus1[2 * l]) /
                                 (pow(2, l) * two_alpha_to_am32));
    return contr.at(c).coeff[p] * one_over_N;
  }

  /// extract primitive shell

  /// @param p the index of the primitive to extract
  /// @param unit_normalized whether to produce unit-normalized primitive; set
  /// to false to produce normalization-free primitive
  /// @return a primitive Shell
  Shell extract_primitive(size_t p, bool unit_normalized = true) const {
    assert(p < nprim());
    svector<Contraction> prim_contr;
    prim_contr.reserve(ncontr());
    for (auto&& c : contr) {
      prim_contr.emplace_back(Contraction{c.l, c.pure, {1.}});
    }
    return Shell({alpha[p]}, prim_contr, O, unit_normalized);
  }

 private:
  // this makes a unit shell
  struct make_unit {};
  Shell(make_unit)
      : alpha{0.0},              // exponent = 0
        contr{{0, false, {1}}},  // contraction coefficient = 1
        O{{0, 0, 0}},            // placed at origin
        max_ln_coeff{0} {}

  /// embeds normalization constants into contraction coefficients. Do this
  /// before computing integrals. \warning Must be done only once. \note this is
  /// now private
  void renorm() {
    using libint2::math::df_Kminus1;
    using std::pow;
    const auto sqrt_Pi_cubed = real_t{5.56832799683170784528481798212};
    const auto np = nprim();
    for (auto& c : contr) {
      assert(c.l <= 15);  // due to df_Kminus1[] a 64-bit integer type; kinda
                          // ridiculous restriction anyway
      for (auto p = 0ul; p != np; ++p) {
        assert(alpha[p] >= 0);
        if (alpha[p] != 0) {
          const auto two_alpha = 2 * alpha[p];
          const auto two_alpha_to_am32 =
              pow(two_alpha, c.l + 1) * sqrt(two_alpha);
          const auto normalization_factor =
              sqrt(pow(2, c.l) * two_alpha_to_am32 /
                   (sqrt_Pi_cubed * df_Kminus1[2 * c.l]));

          c.coeff[p] *= normalization_factor;
        }
      }

      // need to force normalization to unity?
      if (do_enforce_unit_normalization()) {
        // compute the self-overlap of the , scale coefficients by its inverse
        // square root
        double norm{0};
        for (auto p = 0ul; p != np; ++p) {
          for (decltype(p) q = 0ul; q <= p; ++q) {
            auto gamma = alpha[p] + alpha[q];
            norm += (p == q ? 1 : 2) * df_Kminus1[2 * c.l] * sqrt_Pi_cubed *
                    c.coeff[p] * c.coeff[q] /
                    (pow(2, c.l) * pow(gamma, c.l + 1) * sqrt(gamma));
          }
        }
        auto normalization_factor = 1 / sqrt(norm);
        for (auto p = 0ul; p != np; ++p) {
          c.coeff[p] *= normalization_factor;
        }
      }
    }

    update_max_ln_coeff();
  }

  void update_max_ln_coeff() {
    // update max log coefficients
    max_ln_coeff.resize(nprim());
    for (auto p = 0ul; p != nprim(); ++p) {
      real_t max_ln_c = -std::numeric_limits<real_t>::max();
      for (auto& c : contr) {
        max_ln_c = std::max(max_ln_c, std::log(std::abs(c.coeff[p])));
      }
      max_ln_coeff[p] = max_ln_c;
    }
  }
};

inline std::ostream& operator<<(std::ostream& os, const Shell& sh) {
  os << "Shell:( O={" << sh.O[0] << "," << sh.O[1] << "," << sh.O[2] << "}"
     << std::endl;
  os << "  ";
  for (const auto& c : sh.contr) {
    os << " {l=" << c.l << ",sph=" << c.pure << "}";
  }
  os << std::endl;

  for (auto i = 0ul; i < sh.alpha.size(); ++i) {
    os << "  " << sh.alpha[i];
    for (const auto& c : sh.contr) {
      os << " " << c.coeff.at(i);
    }
    os << std::endl;
  }

  return os;
}

// clang-format off
  /// @brief describes method for primitive screening used by ShellPair and Engine
  ///
  /// @note *Rationale*. Since Engine::compute2() can compute primitive data on-the-fly it needs cheap methods for screening primitives.
  ///       - ScreeningMethod::Original was the fast approach that works well for uncontracted integrals over spherical (L=0) Gaussians, but can introduce significant errors in integrals over nonspherical and contracted Gaussians.
  ///       - ScreeningMethod::Conservative addresses the weaknesses of the original method and works well across the board.
  ///       - ScreeningMethod::Schwarz and ScreeningMethod::SchwarzInf should be used when it is possible to precompute the shell pair data _and_ the data can be computed for a specific operator type. These approaches provide rigorous guarantees of precision.
  enum class ScreeningMethod {
    /// standard screening method:
    /// - omit primitive pair if \f$ {\rm ln_scr} \equiv \max|c_a| \max|c_b| \exp(- \alpha_a \alpha_b |AB|^2 / (\alpha_a + \alpha_b) ) < \epsilon \f$
    /// - omit primitive 2-body integral if `spbra.primpairs[pb].ln_scr + spket.primpairs[pk].ln_scr < log(ε)` or `scale * c_a * c_b * c_c * c_d * spbrapp.K * spketpp.K / sqrt(spbrapp.gamma + spketpp.gamma) < ε`
    Original = 0x0001,
    /// conservative screening:
    /// - omit primitive pair if \f$ {\rm ln_scr} \equiv (k_a k_b) \max|c_a| \max|c_b| \exp(- \rho |AB|^2 ) \sqrt{2 \pi^{5/2}} \max( (\max_i|PA_i|)^{L_a} (\max_i|PB_i|)^{L_b}, L_a! L_b! / (\alpha_a + \alpha_b)^{L_a+L_b} ) / (\alpha_a + \alpha_b) < \epsilon \f$
    /// - omit primitive 2-body integral if `spbra.primpairs[pb].ln_scr + spket.primpairs[pk].ln_scr < log(ε)` or `scale * c_a * c_b * c_c * c_d * max(1, spbrapp.nonsph_screen_fac * spketpp.nonsph_screen_fac) * spbrapp.K * spketpp.K / sqrt(spbrapp.gamma + spketpp.gamma) < ε / (npbra * npket)`
    Conservative = 0x0010,
    /// Schwarz screening method using Frobenius norm:
    /// - omit primitive pair if \f$ {\rm ln_scr} \equiv (k_a k_b) \max|c_a| \max|c_b| K_{ab}  < \epsilon \f$, where \f$ K_{ab} \equiv \sqrt{||(ab|ab)||_2} \f$
    /// - omit primitive 2-body integral if `spbra.primpairs[pb].ln_scr + spket.primpairs[pk].ln_scr < log(ε)`
    Schwarz = 0x0100,
    /// Schwarz screening method using infinity norm:
    /// - omit primitive pair if \f$ {\rm ln_scr} \equiv (k_a k_b)  \max|c_a| \max|c_b| K_{ab}  < \epsilon \f$, where \f$ K_{ab} \equiv \sqrt{\max{|(ab|ab)|}} \f$
    /// - omit primitive 2-body integral if `spbra.primpairs[pb].ln_scr + spket.primpairs[pk].ln_scr < log(ε)`
    SchwarzInf = 0x1000,
    Invalid = 0x0000
  };
// clang-format on

namespace detail {
inline ScreeningMethod& default_screening_method_accessor() {
  static ScreeningMethod default_screening_method = ScreeningMethod::Original;
  return default_screening_method;
}
}  // namespace detail

inline ScreeningMethod default_screening_method() {
  return detail::default_screening_method_accessor();
}

inline void default_screening_method(ScreeningMethod screening_method) {
  detail::default_screening_method_accessor() = screening_method;
}

/// ShellPair contains pre-computed shell-pair data, primitive pairs are
/// screened to finite precision
struct ShellPair {
  typedef Shell::real_t real_t;

  // clang-format off
      /// PrimPairData contains pre-computed primitive pair data
      struct PrimPairData {
          real_t P[3]; //!< \f$ (\alpha_1 \vec{A} + \alpha_2 \vec{B})/(\alpha_1 + \alpha_2) \f$
          real_t K;    //!< \f$ \sqrt{2} \pi^{5/4} \exp(-|\vec{A}-\vec{B}|^2 \alpha_a * \alpha_b / (\alpha_a + \alpha_b)) / (\alpha_a + \alpha_b)  \f$
          real_t one_over_gamma; //!< \f$ 1 / (\alpha_a + \alpha_b)  \f$
          real_t nonsph_screen_fac; //!< used only when `screening_method_==ScreeningMethod::Conservative`: approximate upper bound for the modulation of the integrals due to nonspherical bra: \f$  \max( (\max_i|PA_i|)^{L_a} (\max_i|PB_i|)^{L_b}, L_a! L_b! / (\alpha_a + \alpha_b)^{L_a+L_b} ) \f$
          real_t ln_scr; //!< natural log of the primitive pair screening factor (see ScreeningMethod )
          int p1;  //!< first primitive index
          int p2;  //!< second primitive index
      };
  // clang-format on

  std::vector<PrimPairData> primpairs;
  real_t AB[3];
  real_t ln_prec = std::numeric_limits<real_t>::lowest();
  ScreeningMethod screening_method_ = ScreeningMethod::Invalid;

  ShellPair() : primpairs() {
    for (int i = 0; i != 3; ++i) AB[i] = 0.;
  }

  ShellPair(size_t max_nprim) : primpairs() {
    primpairs.reserve(max_nprim * max_nprim);
    for (int i = 0; i != 3; ++i) AB[i] = 0.;
  }
  template <typename Real>
  ShellPair(const Shell& s1, const Shell& s2, Real ln_prec,
            ScreeningMethod screening_method = default_screening_method()) {
    init(s1, s2, ln_prec, screening_method);
  }
  template <typename Real, typename SchwarzFactorEvaluator>
  ShellPair(const Shell& s1, const Shell& s2, Real ln_prec,
            ScreeningMethod screening_method,
            SchwarzFactorEvaluator&& schwarz_factor_evaluator) {
    init(s1, s2, ln_prec, screening_method,
         std::forward<SchwarzFactorEvaluator>(schwarz_factor_evaluator));
  }

  /// makes this equivalent to a default-initialized ShellPair, however the
  /// memory allocated in primpairs is not released
  void reset() {
    primpairs.clear();
    for (int i = 0; i != 3; ++i) AB[i] = 0.;
    ln_prec = std::numeric_limits<real_t>::lowest();
    screening_method_ = ScreeningMethod::Invalid;
  }

  void resize(std::size_t max_nprim) {
    const auto max_nprim2 = max_nprim * max_nprim;
    if (max_nprim * max_nprim > primpairs.size()) primpairs.resize(max_nprim2);
  }

  /// initializes the shell pair data using original or conservative screening
  /// methods
  template <typename Real>
  void init(const Shell& s1, const Shell& s2, Real ln_prec,
            ScreeningMethod screening_method = ScreeningMethod::Original) {
    assert(screening_method == ScreeningMethod::Original ||
           screening_method == ScreeningMethod::Conservative);

    using std::log;

    primpairs.clear();

    const auto& A = s1.O;
    const auto& B = s2.O;
    real_t AB2 = 0.;
    for (int i = 0; i != 3; ++i) {
      AB[i] = A[i] - B[i];
      AB2 += AB[i] * AB[i];
    }

    auto max_l = [](const Shell& s) {
      using std::begin;
      using std::end;
      return std::max_element(
                 begin(s.contr), end(s.contr),
                 [](const Shell::Contraction& c1,
                    const Shell::Contraction& c2) { return c1.l < c2.l; })
          ->l;
    };
    const auto max_l1 = max_l(s1);
    const auto max_l2 = max_l(s2);

    const auto nprim1 = s1.alpha.size();
    const auto nprim2 = s2.alpha.size();
    size_t c = 0;
    for (size_t p1 = 0; p1 != nprim1; ++p1) {
      for (size_t p2 = 0; p2 != nprim2; ++p2) {
        const auto& a1 = s1.alpha[p1];
        const auto& a2 = s2.alpha[p2];
        const auto gamma = a1 + a2;
        const auto oogamma = 1 / gamma;

        const auto rho = a1 * a2 * oogamma;
        const auto minus_rho_times_AB2 = -rho * AB2;
        real_t ln_screen_fac =
            minus_rho_times_AB2 + s1.max_ln_coeff[p1] + s2.max_ln_coeff[p2];
        if (screening_method == ScreeningMethod::Original &&
            ln_screen_fac < ln_prec)
          continue;

        real_t P[3];
        if (AB2 == 0.) {  // this buys a bit more precision
          P[0] = A[0];
          P[1] = A[1];
          P[2] = A[2];
        } else {
          P[0] = (a1 * A[0] + a2 * B[0]) * oogamma;
          P[1] = (a1 * A[1] + a2 * B[1]) * oogamma;
          P[2] = (a1 * A[2] + a2 * B[2]) * oogamma;
        }

        // conservative screening:
        // - partitions the error among all primitive pairs (use \epsilon /
        // nprim to screen, instead of \epsilon itself), and
        // - accounts for the proper spherical gaussian prefactor in the
        // integrals (namely, adds extra \sqrt{2 \pi^{52}}/\gamma_{ab} factor)
        // - accounts for the nonspherical gaussians ... namely
        //   magnitude of primitive (ab|00) integral for nonzero L differs from
        //   that of (00|00) by the magnitude of:
        //   - (max_i|PA_i|)^La (max_i|PB_i|)^Lb when A-B separation is large,
        //   or
        //   - La! Lb! / gammap^(La+Lb) when the separation is small
        real_t nonspherical_scr_factor = 0;
        if (screening_method == ScreeningMethod::Conservative) {
          const auto maxabs_PA_i_to_l1 = std::pow(
              std::max(std::max(std::abs(P[0] - A[0]), std::abs(P[1] - A[1])),
                       std::abs(P[2] - A[2])),
              max_l1);
          const auto maxabs_PB_i_to_l2 = std::pow(
              std::max(std::max(std::abs(P[0] - B[0]), std::abs(P[1] - B[1])),
                       std::abs(P[2] - B[2])),
              max_l2);
          const auto fac_l1_fac_l2_oogamma_to_l =
              math::fac[max_l1] * math::fac[max_l2] *
              std::pow(oogamma, max_l1 + max_l2);
          nonspherical_scr_factor =
              std::max(maxabs_PA_i_to_l1 * maxabs_PB_i_to_l2,
                       fac_l1_fac_l2_oogamma_to_l);
          const auto ln_nonspherical_scr_factor =
              log(std::max(nonspherical_scr_factor, static_cast<real_t>(1)));

          constexpr decltype(rho) ln_sqrt_two_times_M_PI_to_1pt25 =
              1.777485947591722872387900;  // \ln(\sqrt{2} (\pi)^{5/4})
          const auto ln_spherical_scr_extra_factor =
              ln_sqrt_two_times_M_PI_to_1pt25 + log(oogamma);
          const auto ln_nprim = log(nprim1 * nprim2);
          ln_screen_fac += ln_spherical_scr_extra_factor +
                           ln_nonspherical_scr_factor + ln_nprim;
          if (ln_screen_fac < ln_prec) continue;
        }

        primpairs.resize(c + 1);
        PrimPairData& p = primpairs[c];
        p.ln_scr = ln_screen_fac;
        p.p1 = p1;
        p.p2 = p2;
        constexpr decltype(rho) sqrt_two_times_M_PI_to_1pt25 =
            5.9149671727956128778;  // \sqrt{2} (\pi)^{5/4}
        p.K = sqrt_two_times_M_PI_to_1pt25 * exp(minus_rho_times_AB2) * oogamma;
        p.P[0] = P[0];
        p.P[1] = P[1];
        p.P[2] = P[2];
        p.nonsph_screen_fac = nonspherical_scr_factor;
        p.one_over_gamma = oogamma;

        ++c;
      }
    }

    this->ln_prec = ln_prec;
    this->screening_method_ = screening_method;
  }

  /// initializes the shell pair data using Schwarz screening methods
  template <typename Real, typename SchwarzFactorEvaluator>
  void init(const Shell& s1, const Shell& s2, Real ln_prec,
            ScreeningMethod screening_method,
            SchwarzFactorEvaluator&& schwarz_factor_evaluator) {
    assert(screening_method == ScreeningMethod::Schwarz ||
           screening_method == ScreeningMethod::SchwarzInf);

    using std::log;

    primpairs.clear();

    const auto& A = s1.O;
    const auto& B = s2.O;
    real_t AB2 = 0.;
    for (int i = 0; i != 3; ++i) {
      AB[i] = A[i] - B[i];
      AB2 += AB[i] * AB[i];
    }

    const auto nprim1 = s1.alpha.size();
    const auto nprim2 = s2.alpha.size();
    const auto nprim12 = nprim1 * nprim2;
    size_t c = 0;
    for (size_t p1 = 0; p1 != nprim1; ++p1) {
      for (size_t p2 = 0; p2 != nprim2; ++p2) {
        const auto ln_screen_fac =
            log(nprim12 * schwarz_factor_evaluator(s1, p1, s2, p2)) +
            s1.max_ln_coeff[p1] + s2.max_ln_coeff[p2];
        if (ln_screen_fac < ln_prec) continue;

        const auto& a1 = s1.alpha[p1];
        const auto& a2 = s2.alpha[p2];
        const auto gamma = a1 + a2;
        const auto oogamma = 1 / gamma;

        const auto rho = a1 * a2 * oogamma;
        const auto minus_rho_times_AB2 = -rho * AB2;

        real_t P[3];
        if (AB2 == 0.) {  // this buys a bit more precision
          P[0] = A[0];
          P[1] = A[1];
          P[2] = A[2];
        } else {
          P[0] = (a1 * A[0] + a2 * B[0]) * oogamma;
          P[1] = (a1 * A[1] + a2 * B[1]) * oogamma;
          P[2] = (a1 * A[2] + a2 * B[2]) * oogamma;
        }

        primpairs.resize(c + 1);
        PrimPairData& p = primpairs[c];
        p.ln_scr = ln_screen_fac;
        p.p1 = p1;
        p.p2 = p2;
        constexpr decltype(rho) sqrt_two_times_M_PI_to_1pt25 =
            5.9149671727956128778;  // \sqrt{2} (\pi)^{5/4}
        p.K = sqrt_two_times_M_PI_to_1pt25 * exp(minus_rho_times_AB2) * oogamma;
        p.P[0] = P[0];
        p.P[1] = P[1];
        p.P[2] = P[2];
        p.nonsph_screen_fac = 0;
        p.one_over_gamma = oogamma;

        ++c;
      }
    }

    this->ln_prec = ln_prec;
    this->screening_method_ = screening_method;
  }
};

}  // namespace libint2

#endif /* _libint2_src_lib_libint_shell_h_ */
