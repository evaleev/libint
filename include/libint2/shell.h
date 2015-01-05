/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License, version 2,
 *  as published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

#ifndef _libint2_src_lib_libint_shell_h_
#define _libint2_src_lib_libint_shell_h_

#if __cplusplus <= 199711L
# error "The simple Libint API requires C++11 support"
#endif

#include <iostream>
#include <array>
#include <vector>

#include <libint2.h>

namespace libint2 {

  namespace math {
    /// fac[k] = k!
    static constexpr std::array<int64_t,21> fac = {{1L, 1L, 2L, 6L, 24L, 120L, 720L, 5040L, 40320L, 362880L, 3628800L, 39916800L,
                                                    479001600L, 6227020800L, 87178291200L, 1307674368000L, 20922789888000L,
                                                    355687428096000L, 6402373705728000L, 121645100408832000L,
                                                    2432902008176640000L}};
    /// df_Kminus1[k] = (k-1)!!
    static constexpr std::array<int64_t,31> df_Kminus1 = {{1, 1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080, 135135,
                                                           645120, 2027025, 10321920, 34459425, 185794560, 654729075,
                                                           3715891200, 13749310575, 81749606400, 316234143225, 1961990553600,
                                                           7905853580625, 51011754393600, 213458046676875, 1428329123020800,
                                                           6190283353629375}};
    /// bc(i,j) = binomial coefficient, i! / (j! (i-j)!)
    template <typename Int> int64_t bc(Int i, Int j) {
      assert(i < fac.size());
      assert(j < fac.size());
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
   *                    {0, false, {-0.09996723, 0.39951283, 0.70011547}}
   *                    {1, false, {0.15591627, 0.60768372, 0.39195739}}
   *                  },
   *                   {{0.0, 0.0, 0.0}}
   *                };
   *  \endverbatim
   *  \note The contraction coefficients correspond to <em>unity-normalized</em> primitives. EMSL Gaussian Basis Set Database,
   *  as well as basis set libraries embedded into most quantum chemistry programs use this convention.
   *  However, coefficients must be converted to refer to normalization-free primitives before computing
   *  integrals! See Shell::renorm();
   */
  struct Shell {
      typedef ::libint2::real_t real_t;
      typedef ::libint2::realvec_t realvec_t;

      /// contracted Gaussian = angular momentum + sph/cart flag + contraction coefficients
      struct Contraction {
          int l;
          bool pure;
          std::vector<real_t> coeff;
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

      std::vector<real_t> alpha; //!< exponents
      std::vector<Contraction> contr;      //!< contractions
      std::array<real_t, 3> O;   //!< origin
      std::vector<real_t> max_ln_coeff; //!< maximum ln of (absolute) contraction coefficient for each primitive

      Shell& move(const std::array<real_t, 3> new_origin) {
        O = new_origin;
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

      /// embeds normalization constants into contraction coefficients. Do this before computing integrals.
      /// \note Must be done only once.
      void renorm() {
        using libint2::math::df_Kminus1;
        const auto sqrt_Pi_cubed = real_t{5.56832799683170784528481798212};
        const auto np = nprim();
        for(auto& c: contr) {
          assert(c.l <= 15); // due to df_Kminus1[] a 64-bit integer type; kinda ridiculous restriction anyway
          for(auto p=0; p!=np; ++p) {
            assert(alpha[p] >= 0.0);
            if (alpha[p] != 0.) {
              const auto two_alpha = 2 * alpha[p];
              const auto two_alpha_to_am32 = pow(two_alpha,c.l+1) * sqrt(two_alpha);
              const auto norm = sqrt(pow(2,c.l) * two_alpha_to_am32/(sqrt_Pi_cubed * df_Kminus1[2*c.l] ));

              c.coeff[p] *= norm;
            }
          }
        }

        // update max log coefficients
        max_ln_coeff.resize(np);
        for(auto p=0; p!=np; ++p) {
          real_t max_ln_c = - std::numeric_limits<real_t>::max();
          for(auto& c: contr) {
            max_ln_c = std::max(max_ln_c, log(std::abs(c.coeff[p])));
          }
          max_ln_coeff[p] = max_ln_c;
        }
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
  };

  inline std::ostream& operator<<(std::ostream& os, const Shell& sh) {
    os << "Shell:( O={" << sh.O[0] << "," << sh.O[1] << "," << sh.O[2] << "}" << std::endl;
    os << "  ";
    for(const auto& c: sh.contr) {
      os << " {l=" << c.l << ",sph=" << c.pure << "}";
    }
    os << std::endl;

    for(auto i=0; i<sh.alpha.size(); ++i) {
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

      // initializes "expensive" primitive pair data; primitive pairs are
      void init(const Shell& s1, const Shell& s2, const real_t& ln_prec) {

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
            const auto oogamma = 1.0 / gamma;

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
            p.P[0] = (a1 * A[0] + a2 * B[0]) * oogamma;
            p.P[1] = (a1 * A[1] + a2 * B[1]) * oogamma;
            p.P[2] = (a1 * A[2] + a2 * B[2]) * oogamma;
            p.one_over_gamma = oogamma;

            ++c;
          }
        }
      }

  };

} // namespace libint2

#endif /* _libint2_src_lib_libint_shell_h_ */
