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

  /** Simple-to-use Gaussian shell. Here's an example of how to create an s+p shell of the STO-3G basis on the oxygen atom
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
   *  \note The contraction coefficients correspond to <em>unity-normalized</em> primitives. Basis Set Database,
   *  as well as basis set libraries embedded into most quantum chemistry programs use this convention.
   *  However, coefficients must be converted to refer to normalization-free primitives before computing
   *  integrals! See Shell::renorm();
   */
  struct Shell {
      /// contracted Gaussian = angular momentum + sph/cart flag + contraction coefficients
      struct Contraction {
          int l;
          bool pure;
          std::vector<LIBINT2_REALTYPE> coeff;
      };

      std::vector<LIBINT2_REALTYPE> alpha; //!< exponents
      std::vector<Contraction> contr;      //!< contractions
      std::array<LIBINT2_REALTYPE, 3> O;   //!< origin

      Shell& move(const std::array<LIBINT2_REALTYPE, 3> new_origin) {
        O = new_origin;
        return *this;
      }

      size_t size() const {
        size_t s = 0;
        for(const auto& c: contr) {
          s += c.pure ? (2 * c.l + 1) : ((c.l + 1) * (c.l + 2) / 2);
        }
        return s;
      }

      /// embeds normalization constants into contraction coefficients. Do this before computing integrals.
      /// \note Must be done only once.
      void renorm() {
        using libint2::math::df_Kminus1;
        const auto sqrt_Pi_cubed = LIBINT2_REALTYPE{5.56832799683170784528481798212};
        const auto np = nprim();
        for(auto& c: contr) {
          assert(c.l <= 15); // due to df_Kminus1[] a 64-bit integer type; kinda ridiculous restriction anyway
          for(auto p=0; p!=np; ++p) {
            const auto two_alpha = 2 * alpha[p];
            const auto two_alpha_to_am32 = pow(two_alpha,c.l+1) * sqrt(two_alpha);
            const auto norm = sqrt(pow(2,c.l) * two_alpha_to_am32/(sqrt_Pi_cubed * df_Kminus1[2*c.l] ));

            c.coeff[p] *= norm;

          }
        }
      }

      size_t ncontr() const { return contr.size(); }
      size_t nprim() const { return alpha.size(); }

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

} // namespace libint2

#endif /* _libint2_src_lib_libint_shell_h_ */
