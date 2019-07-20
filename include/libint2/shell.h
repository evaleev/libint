/*
 *  Copyright (C) 2004-2019 Edward F. Valeev
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
# error "libint2/shell.h requires C++11 support"
#endif

#include <iostream>
#include <array>
#include <vector>
#include <cassert>
#include <cmath>

#include <libint2.h>

#include <libint2/util/small_vector.h>

namespace libint2 {

  namespace math {
    /// fac[k] = k!
    static constexpr std::array<int64_t,21> fac = {{1LL, 1LL, 2LL, 6LL, 24LL, 120LL, 720LL, 5040LL, 40320LL, 362880LL, 3628800LL, 39916800LL,
                                                    479001600LL, 6227020800LL, 87178291200LL, 1307674368000LL, 20922789888000LL,
                                                    355687428096000LL, 6402373705728000LL, 121645100408832000LL,
                                                    2432902008176640000LL}};
    /// df_Kminus1[k] = (k-1)!!
    static constexpr std::array<int64_t,31> df_Kminus1 = {{1LL, 1LL, 1LL, 2LL, 3LL, 8LL, 15LL, 48LL, 105LL, 384LL, 945LL, 3840LL, 10395LL, 46080LL, 135135LL,
                                                           645120LL, 2027025LL, 10321920LL, 34459425LL, 185794560LL, 654729075LL,
                                                           3715891200LL, 13749310575LL, 81749606400LL, 316234143225LL, 1961990553600LL,
                                                           7905853580625LL, 51011754393600LL, 213458046676875LL, 1428329123020800LL,
                                                           6190283353629375LL}};
    /// bc(i,j) = binomial coefficient, i! / (j! (i-j)!)
    template <typename Int> int64_t bc(Int i, Int j) {
      assert(i < Int(fac.size()));
      assert(j < Int(fac.size()));
      assert(i >= j);
      return fac[i] / (fac[j] * fac[i-j]);
    }
  }

  /// generally-contracted Solid-Harmonic/Cartesion Gaussian Shell

  /** A simple-to-use Gaussian shell. Here's an example of how to create an s+p shell of the STO-3G basis on the oxygen atom
   *  located at the origin.
   *  \verbatim
   *  auto s = Shell{
   *                  {5.033151300, 1.169596100, 0.380389000},
   *                  {
   *                    {0, false, {-0.09996723, 0.39951283, 0.70011547}},
   *                    {1, false, {0.15591627, 0.60768372, 0.39195739}}
   *                  },
   *                   {{0.0, 0.0, 0.0}}
   *                };
   *  \endverbatim
   *  \note The contraction coefficients correspond to <em>unity-normalized</em> primitives. EMSL Gaussian Basis Set Database,
   *  as well as basis set libraries embedded into most quantum chemistry programs use this convention.
   *  However, coefficients are automatically converted internally to refer to
   *  normalization-free primitives before computing integrals (see \c Shell::renorm() ).
   */
  struct Shell {
      typedef double real_t;

      /// contracted Gaussian = angular momentum + sph/cart flag + contraction coefficients
      struct Contraction {
          int l;
          bool pure;
          svector<real_t> coeff;

          bool operator==(const Contraction& other) const {
            return &other == this || (l == other.l && pure == other.pure && coeff == other.coeff);
          }
          bool operator!=(const Contraction& other) const {
            return not this->operator==(other);
          }
          size_t cartesian_size() const {
            return (l + 1) * (l + 2) / 2;
          }
          size_t size() const {
              return pure ? (2 * l + 1) : cartesian_size();
          }
      };

      svector<real_t> alpha; //!< exponents
      svector<Contraction> contr;      //!< contractions
      std::array<real_t, 3> O;   //!< origin
      svector<real_t> max_ln_coeff; //!< maximum ln of (absolute) contraction coefficient for each primitive

      Shell() = default;
      Shell(const Shell&) = default;
      // intel does not support "move ctor = default"
      Shell(Shell&& other) noexcept :
        alpha(std::move(other.alpha)),
        contr(std::move(other.contr)),
        O(std::move(other.O)),
        max_ln_coeff(std::move(other.max_ln_coeff)) {
      }
      Shell& operator=(const Shell&) = default;
      // intel does not support "move asgnmt = default"
      Shell& operator=(Shell&& other) noexcept {
        alpha = std::move(other.alpha);
        contr = std::move(other.contr);
        O = std::move(other.O);
        max_ln_coeff = std::move(other.max_ln_coeff);
        return *this;
      }
      Shell(svector<real_t> _alpha,
            svector<Contraction> _contr,
            std::array<real_t, 3> _O) :
              alpha(std::move(_alpha)),
              contr(std::move(_contr)),
              O(std::move(_O)) {
        // embed normalization factors into contraction coefficients
        renorm();
      }

      Shell& move(std::array<real_t, 3> new_origin) {
        O = std::move(new_origin);
        return *this;
      }

      size_t cartesian_size() const {
        size_t s = 0;
        for(const auto& c: contr) { s += c.cartesian_size(); }
        return s;
      }
      size_t size() const {
        size_t s = 0;
        for(const auto& c: contr) { s += c.size(); }
        return s;
      }

      size_t ncontr() const { return contr.size(); }
      size_t nprim() const { return alpha.size(); }

      bool operator==(const Shell& other) const {
        return &other == this || (O == other.O && alpha == other.alpha && contr == other.contr);
      }
      bool operator!=(const Shell& other) const {
        return not this->operator==(other);
      }

      static char am_symbol(size_t l) {
        static char lsymb[] = "spdfghikmnoqrtuvwxyz";
        assert(l<=19);
        return lsymb[l];
      }
      static unsigned short am_symbol_to_l(char am_symbol) {
        const char AM_SYMBOL = ::toupper(am_symbol);
        switch (AM_SYMBOL) {
          case 'S': return 0;
          case 'P': return 1;
          case 'D': return 2;
          case 'F': return 3;
          case 'G': return 4;
          case 'H': return 5;
          case 'I': return 6;
          case 'K': return 7;
          case 'M': return 8;
          case 'N': return 9;
          case 'O': return 10;
          case 'Q': return 11;
          case 'R': return 12;
          case 'T': return 13;
          case 'U': return 14;
          case 'V': return 15;
          case 'W': return 16;
          case 'X': return 17;
          case 'Y': return 18;
          case 'Z': return 19;
          default: throw "invalid angular momentum label";
        }
      }

      struct defaultable_boolean {
          typedef enum {false_value=0,true_value=1,default_value=2} value_t;
          defaultable_boolean() : value_(default_value) {}
          defaultable_boolean(bool v) : value_(static_cast<value_t>(v?1:0)) {}
          bool is_default() const { return value_ == default_value; }
          operator bool() const { assert(value_ != default_value); return value_ == true_value; }
        private:
          value_t value_;
      };

      /// sets and/or reports whether the auto-renormalization to unity is set
      /// if called without arguments, returns the current value of the flag
      /// otherwise, will set the flag to \c flag
      /// \note by default, shells WILL be re-normalized to unity
      static bool do_enforce_unit_normalization(defaultable_boolean flag = defaultable_boolean()) {
        static bool result{true};
        if (not flag.is_default()) {
          result = flag;
        }
        return result;
      }

      /// @return "unit" Shell, with exponent=0. and coefficient=1., located at the origin
      static const Shell& unit() {
        static const Shell unitshell{make_unit()};
        return unitshell;
      }

      /// @return the coefficient of primitive \c p in contraction \c c assuming unit normalized primitive
      ///         (coeff contains coefficients of normalization-free primitives; @sa Shell::renorm() )
      real_t coeff_normalized(size_t c, size_t p) const {
        const auto alpha = this->alpha.at(p);
        assert(alpha >= 0.0);
        const auto l = contr.at(c).l;
        assert(l <= 15); // due to df_Kminus1[] a 64-bit integer type; kinda ridiculous restriction anyway

        using libint2::math::df_Kminus1;
        const auto sqrt_Pi_cubed = real_t{5.56832799683170784528481798212};
        const auto two_alpha = 2 * alpha;
        const auto two_alpha_to_am32 = pow(two_alpha,l+1) * sqrt(two_alpha);
        const auto one_over_N = sqrt((sqrt_Pi_cubed * df_Kminus1[2*l] )/(pow(2,l) * two_alpha_to_am32));
        return contr.at(c).coeff[p] * one_over_N;
      }

    private:

      // this makes a unit shell
      struct make_unit{};
      Shell(make_unit) :
        alpha{0.0},                           // exponent = 0
        contr{{0, false, {1}}},  // contraction coefficient = 1
        O{{0, 0, 0}},                   // placed at origin
        max_ln_coeff{0} {
      }

      /// embeds normalization constants into contraction coefficients. Do this before computing integrals.
      /// \warning Must be done only once.
      /// \note this is now private
      void renorm() {
        using libint2::math::df_Kminus1;
        using std::pow;
        const auto sqrt_Pi_cubed = real_t{5.56832799683170784528481798212};
        const auto np = nprim();
        for(auto& c: contr) {
          assert(c.l <= 15); // due to df_Kminus1[] a 64-bit integer type; kinda ridiculous restriction anyway
          for(auto p=0ul; p!=np; ++p) {
            assert(alpha[p] >= 0);
            if (alpha[p] != 0) {
              const auto two_alpha = 2 * alpha[p];
              const auto two_alpha_to_am32 = pow(two_alpha,c.l+1) * sqrt(two_alpha);
              const auto normalization_factor = sqrt(pow(2,c.l) * two_alpha_to_am32/(sqrt_Pi_cubed * df_Kminus1[2*c.l] ));

              c.coeff[p] *= normalization_factor;
            }
          }

          // need to force normalization to unity?
          if (do_enforce_unit_normalization()) {
            // compute the self-overlap of the , scale coefficients by its inverse square root
            double norm{0};
            for(auto p=0ul; p!=np; ++p) {
              for(decltype(p) q=0ul; q<=p; ++q) {
                auto gamma = alpha[p] + alpha[q];
                norm += (p==q ? 1 : 2) * df_Kminus1[2*c.l] * sqrt_Pi_cubed * c.coeff[p] * c.coeff[q] /
                        (pow(2,c.l) * pow(gamma,c.l+1) * sqrt(gamma));
              }
            }
            auto normalization_factor = 1 / sqrt(norm);
            for(auto p=0ul; p!=np; ++p) {
              c.coeff[p] *= normalization_factor;
            }
          }

        }

        // update max log coefficients
        max_ln_coeff.resize(np);
        for(auto p=0ul; p!=np; ++p) {
          real_t max_ln_c = - std::numeric_limits<real_t>::max();
          for(auto& c: contr) {
            max_ln_c = std::max(max_ln_c, std::log(std::abs(c.coeff[p])));
          }
          max_ln_coeff[p] = max_ln_c;
        }
      }

  };

  inline std::ostream& operator<<(std::ostream& os, const Shell& sh) {
    os << "Shell:( O={" << sh.O[0] << "," << sh.O[1] << "," << sh.O[2] << "}" << std::endl;
    os << "  ";
    for(const auto& c: sh.contr) {
      os << " {l=" << c.l << ",sph=" << c.pure << "}";
    }
    os << std::endl;

    for(auto i=0ul; i<sh.alpha.size(); ++i) {
      os << "  " << sh.alpha[i];
      for(const auto& c: sh.contr) {
        os << " " << c.coeff.at(i);
      }
      os << std::endl;
    }

    return os;
  }

  /// ShellPair pre-computes shell-pair data, primitive pairs are screened to finite precision
  struct ShellPair {
      typedef Shell::real_t real_t;

      struct PrimPairData {
          real_t P[3]; //!< \f$ (\alpha_1 \vec{A} + \alpha_2 \vec{B})/(\alpha_1 + \alpha_2) \f$
          real_t K;
          real_t one_over_gamma;
          real_t scr;
          int p1;
          int p2;
      };

      std::vector<PrimPairData> primpairs;
      real_t AB[3];

      ShellPair() : primpairs() { for(int i=0; i!=3; ++i) AB[i] = 0.; }

      ShellPair(size_t max_nprim) : primpairs() {
        primpairs.reserve(max_nprim*max_nprim);
        for(int i=0; i!=3; ++i) AB[i] = 0.;
      }
      template <typename Real> ShellPair(const Shell& s1, const Shell& s2, Real ln_prec) {
        init(s1, s2, ln_prec);
      }

      void resize(std::size_t max_nprim) {
        const auto max_nprim2 = max_nprim * max_nprim;
        if (max_nprim * max_nprim > primpairs.size())
          primpairs.resize(max_nprim2);
      }

      /// initializes "expensive" primitive pair data; a pair of primitives with exponents \f$ \{\alpha_a,\alpha_b\} \f$
      /// located at \f$ \{ \vec{A},\vec{B} \} \f$ whose max coefficients in contractions are \f$ \{ \max{|c_a|} , \max{|c_b|} \} \f$ is screened-out (omitted)
      /// if \f$ \exp(-|\vec{A}-\vec{B}|^2 \alpha_a * \alpha_b / (\alpha_a + \alpha_b)) \max{|c_a|} \max{|c_b|} \leq \epsilon \f$
      /// where \f$ \epsilon \f$ is the desired precision of the integrals.
      template <typename Real> void init(const Shell& s1, const Shell& s2, Real ln_prec) {

        primpairs.clear();

        const auto& A = s1.O;
        const auto& B = s2.O;
        real_t AB2 = 0.;
        for(int i=0; i!=3; ++i) {
          AB[i] = A[i] - B[i];
          AB2 += AB[i]*AB[i];
        }

        size_t c = 0;
        for(size_t p1=0; p1!=s1.alpha.size(); ++p1) {
          for(size_t p2=0; p2!=s2.alpha.size(); ++p2) {

            const auto& a1 = s1.alpha[p1];
            const auto& a2 = s2.alpha[p2];
            const auto gamma = a1 + a2;
            const auto oogamma = 1 / gamma;

            const auto rho = a1 * a2 * oogamma;
            const auto minus_rho_times_AB2 = -rho*AB2;
            const auto screen_fac = minus_rho_times_AB2 + s1.max_ln_coeff[p1] + s2.max_ln_coeff[p2];
            if (screen_fac < ln_prec)
              continue;

            primpairs.resize(c+1);
            PrimPairData& p = primpairs[c];
            p.scr = screen_fac;
            p.p1 = p1;
            p.p2 = p2;
            p.K = exp(minus_rho_times_AB2) * oogamma;
            if (AB2 == 0.) {  // this buys a bit more precision
              p.P[0] = A[0];
              p.P[1] = A[1];
              p.P[2] = A[2];
            } else {
              p.P[0] = (a1 * A[0] + a2 * B[0]) * oogamma;
              p.P[1] = (a1 * A[1] + a2 * B[1]) * oogamma;
              p.P[2] = (a1 * A[2] + a2 * B[2]) * oogamma;
            }
            p.one_over_gamma = oogamma;

            ++c;
          }
        }
      }

  };

} // namespace libint2

#endif /* _libint2_src_lib_libint_shell_h_ */
