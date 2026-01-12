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

#ifndef _libint2_src_lib_libint_vrrgtg1dxxxx_h_
#define _libint2_src_lib_libint_vrrgtg1dxxxx_h_
// #include "../../../../SIMD_wrapped_vector/template/simd_wrapped_vector.hpp"
#include <libint2.h>
#include <util_types.h>

#include <cassert>
#include <cstdlib>

namespace libint2 {

/** builds (ab| GTG_1d |cd), the shell set of 2-dimensional integrals needed for
 *Rys quadrature evaluation of 2-body ints.
 *
 *  @tparam CartesianAxis specifies the Cartesian axis with of the integral,
 *valid values are 0 (x), 1 (y), or 2 (z)
 *  @param src0 (00|GTG_1d|00)
 *
 * @note build (a b|c d) as follows: (1) (a+b 0|c+d 0), (2) (a b|c+d 0), (3) (a
 *b|c d).
 **/
template <unsigned int CartesianAxis, int La, int Lb, int Lc, int Ld,
          bool vectorize>
struct VRR_GTG_1d_xx_xx {
  static void compute(const Libint_t* inteval, VectorSIMD<double, npts>* target,
                      VectorSIMD<double, npts>* src0) {
    enum XYZ { x = 0, y = 1, z = 2 };
    assert(CartesianAxis == x || CartesianAxis == y || CartesianAxis == z);
    // assert(vectorize == false);

    const unsigned int veclen = vectorize ? inteval->veclen : 1;

    // corner case: (00|00)
    if (La == 0 && Lb == 0 && Lc == 0 && Ld == 0) {
      for (unsigned int v = 0; v != veclen; ++v) target[v] = src0[v];
      return;
    }

    //---------------------------------------------
    // Part (1): build (a+b 0|c+d 0)
    //---------------------------------------------

    VectorSIMD<double, npts> apb_0_GTG_cpd_0[La + Lb + 1][Lc + Ld + 1];
    apb_0_GTG_cpd_0[0][0] = src0[0];

    const VectorSIMD<double, npts>*pfac0_0, *pfac0_1;
    const VectorSIMD<double, npts>* pfac1_0 = inteval->R12kG12_pfac1_0;
    const VectorSIMD<double, npts>* pfac1_1 = inteval->R12kG12_pfac1_1;
    const VectorSIMD<double, npts>* pfac2 = inteval->R12kG12_pfac2;
    switch (CartesianAxis) {
      case x:
        pfac0_0 = inteval->R12kG12_pfac0_0_x;
        pfac0_1 = inteval->R12kG12_pfac0_1_x;
        break;
      case y:
        pfac0_0 = inteval->R12kG12_pfac0_0_y;
        pfac0_1 = inteval->R12kG12_pfac0_1_y;
        break;
      case z:
        pfac0_0 = inteval->R12kG12_pfac0_0_z;
        pfac0_1 = inteval->R12kG12_pfac0_1_z;
        break;
      default:
        assert(false);
    }

    // build (0 0|1 0)
    if (Lc + Ld > 0) {
      apb_0_GTG_cpd_0[0][1] = pfac0_1[0] * apb_0_GTG_cpd_0[0][0];
#if LIBiINT2_FLOP_COUNT
      inteval->nflops[0] += 1;
#endif
    }

    // build (0 0|c+d 0)
    if (Lc + Ld > 1) {
      for (int c_plus_d = 1; c_plus_d != Lc + Ld; ++c_plus_d) {
        apb_0_GTG_cpd_0[0][c_plus_d + 1] =
            pfac0_1[0] * apb_0_GTG_cpd_0[0][c_plus_d] +
            c_plus_d * pfac1_1[0] * apb_0_GTG_cpd_0[0][c_plus_d - 1];
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 4 * (Lc + Ld - 1);
#endif
    }

    // build (1 0|0 0)
    if (La + Lb > 0) {
      apb_0_GTG_cpd_0[1][0] = pfac0_0[0] * apb_0_GTG_cpd_0[0][0];
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 1;
#endif
    }

    // build (a+b 0|0 0)
    if (La + Lb > 1) {
      for (int a_plus_b = 1; a_plus_b != La + Lb; ++a_plus_b) {
        apb_0_GTG_cpd_0[a_plus_b + 1][0] =
            pfac0_0[0] * apb_0_GTG_cpd_0[a_plus_b][0] +
            a_plus_b * pfac1_0[0] * apb_0_GTG_cpd_0[a_plus_b - 1][0];
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 4 * (La + Lb - 1);
#endif
    }

    // build (1 0|c+d 0)
    if (La + Lb > 0 && Lc + Ld > 0) {
      for (int c_plus_d = 1; c_plus_d <= Lc + Ld; ++c_plus_d) {
        apb_0_GTG_cpd_0[1][c_plus_d] =
            pfac0_0[0] * apb_0_GTG_cpd_0[0][c_plus_d] +
            c_plus_d * pfac2[0] * apb_0_GTG_cpd_0[0][c_plus_d - 1];
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 4 * (Lc + Ld - 1);
#endif
    }

    // build (a+b 0|c+d 0)
    if (La + Lb > 1 && Lc + Ld > 0) {
      for (int a_plus_b = 1; a_plus_b != La + Lb; ++a_plus_b) {
        for (int c_plus_d = 1; c_plus_d <= Lc + Ld; ++c_plus_d) {
          apb_0_GTG_cpd_0[a_plus_b + 1][c_plus_d] =
              pfac0_0[0] * apb_0_GTG_cpd_0[a_plus_b][c_plus_d] +
              a_plus_b * pfac1_0[0] * apb_0_GTG_cpd_0[a_plus_b - 1][c_plus_d] +
              c_plus_d * pfac2[0] * apb_0_GTG_cpd_0[a_plus_b][c_plus_d - 1];
        }
      }
#if LIBINT2_FLOP_COUNT
      inteval->nflops[0] += 7 * (La + Lb - 1) * (Lc + Ld - 1);
#endif
    }

    //---------------------------------------------
    // Part (2): build (a b|c+d 0)
    //---------------------------------------------

    double AB[1];
    switch (CartesianAxis) {
      case x:
        std::cout << "printing before segfault" << std::endl;
        AB[0] = inteval->AB_x[0];
        break;
      case y:
        AB[0] = inteval->AB_y[0];
        break;
      case z:
        AB[0] = inteval->AB_z[0];
        break;
      default:
        assert(false);
    }

    VectorSIMD<double, npts> a_b_GTG_cpd_0[La + 1][Lb + 1][Lc + Ld + 1];
    for (int c_plus_d = 0; c_plus_d <= Lc + Ld; ++c_plus_d) {
      // copy (a+b 0| to a local 0,a+b buffer
      VectorSIMD<double, npts> b_a_GTG[La + Lb + 1][La + Lb + 1];
      for (int a_plus_b = 0; a_plus_b <= La + Lb; ++a_plus_b) {
        b_a_GTG[0][a_plus_b] = apb_0_GTG_cpd_0[a_plus_b][c_plus_d];
      }
      // use HRR to compute b,a
      for (int b = 1; b <= Lb; ++b) {
        for (int a = 0; a <= La + Lb - b; ++a) {
          b_a_GTG[b][a] = b_a_GTG[b - 1][a + 1] + AB[0] * b_a_GTG[b - 1][a];
        }
#if LIBINT2_FLOP_COUNT
        inteval->nflops[0] += 2 * (La + Lb - b + 1);
#endif
      }
      // copy b,a to (a b|
      for (int b = 0; b <= Lb; ++b) {
        for (int a = 0; a <= La; ++a) {
          a_b_GTG_cpd_0[a][b][c_plus_d] = b_a_GTG[b][a];
        }
      }
    }

    //---------------------------------------------
    // Part (3): build (a b|c d)
    //---------------------------------------------

    double CD[1];
    switch (CartesianAxis) {
      case x:
        CD[0] = inteval->CD_x[0];
        break;
      case y:
        CD[0] = inteval->CD_y[0];
        break;
      case z:
        CD[0] = inteval->CD_z[0];
        break;
      default:
        assert(false);
    }

    VectorSIMD<double, npts>* target_a_b_blk_ptr = target;
    const int Nd = (Ld + 1);
    const int Ncd = (Lc + 1) * Nd;
    for (int a = 0; a <= La; ++a) {
      for (int b = 0; b <= Lb; ++b, target_a_b_blk_ptr += Ncd) {
        // copy |c+d 0) to a local 0,c+d buffer
        VectorSIMD<double, npts> d_c_GTG[Lc + Ld + 1][Lc + Ld + 1];
        for (int c_plus_d = 0; c_plus_d <= Lc + Ld; ++c_plus_d) {
          d_c_GTG[0][c_plus_d] = a_b_GTG_cpd_0[a][b][c_plus_d];
        }
        // use HRR to compute d,c
        for (int d = 1; d <= Ld; ++d) {
          for (int c = 0; c <= Lc + Ld - d; ++c) {
            d_c_GTG[d][c] = d_c_GTG[d - 1][c + 1] + CD[0] * d_c_GTG[d - 1][c];
          }
#if LIBINT2_FLOP_COUNT
          inteval->nflops[0] += 2 * (Lc + Ld - d + 1);
#endif
        }
        // copy d,c to |c d)
        for (int d = 0; d <= Ld; ++d) {
          for (int c = 0, cd = d; c <= Lc; ++c, cd += Nd) {
            target_a_b_blk_ptr[cd] = d_c_GTG[d][c];
          }
        }
      }
    }

    // done
  }
};

};  // namespace libint2

#endif  // header guard
