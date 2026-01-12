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

#include <prefactors.h>
#include <smart_ptr.h>

#include <cstdio>

using namespace libint2;

Prefactors::Prefactors()
    : rho(new rdouble("rho")), one_o_2alphasum(new rdouble("oo2ze")) {
  typedef std::shared_ptr<rdouble> rdptr;

  char XY[np][2] = {"P", "Q"};
  char X[np][2][2] = {{"A", "B"}, {"C", "D"}};

  for (unsigned int p = 0; p < np; p++) {
    for (int braket = 0; braket < 2; braket++) {
      char XY_X_str[20];
      snprintf(XY_X_str, 20, "%s%s", XY[p], X[p][braket]);
      rdptr vXY_X_ptr(new rdouble(XY_X_str));
      vXY_X[p][braket] = vXY_X_ptr;

      char W_XY_str[20];
      snprintf(W_XY_str, 20, "W%s", XY[p]);
      rdptr vW_XY_ptr(new rdouble(W_XY_str));
      vW_XY[p] = vW_XY_ptr;

      char xyz_str[] = "xyz";
      for (int xyz = 0; xyz < 3; xyz++) {
        char XY_X_i_str[20];
        snprintf(XY_X_i_str, 20, "%s_%c", XY_X_str, xyz_str[xyz]);
        rdptr XY_X_i_ptr(new rdouble(XY_X_i_str));
        XY_X[p][braket][xyz] = XY_X_i_ptr;

        char W_XY_i_str[20];
        snprintf(W_XY_i_str, 20, "%s_%c", W_XY_str, xyz_str[xyz]);
        rdptr W_XY_i_ptr(new rdouble(W_XY_i_str));
        W_XY[p][xyz] = W_XY_i_ptr;
      }
    }

    char X_Y_str[20];
    snprintf(X_Y_str, 20, "%s%s", X[p][0], X[p][1]);
    rdptr vX_Y_ptr(new rdouble(X_Y_str));
    vX_Y[p] = vX_Y_ptr;
    const char xyz_str[] = "xyz";
    for (int xyz = 0; xyz < 3; xyz++) {
      char X_Y_i_str[20];
      snprintf(X_Y_i_str, 20, "%s_%c", X_Y_str, xyz_str[xyz]);
      rdptr X_Y_i_ptr(new rdouble(X_Y_i_str));
      X_Y[p][xyz] = X_Y_i_ptr;
    }

    char Y_X_str[20];
    snprintf(Y_X_str, 20, "%s%s", X[p][1], X[p][0]);
    rdptr vY_X_ptr(new rdouble(Y_X_str));
    vY_X[p] = vY_X_ptr;
    for (int xyz = 0; xyz < 3; xyz++) {
      char Y_X_i_str[20];
      snprintf(Y_X_i_str, 20, "%s_%c", Y_X_str, xyz_str[xyz]);
      rdptr Y_X_i_ptr(new rdouble(Y_X_i_str));
      Y_X[p][xyz] = Y_X_i_ptr;
    }
  }

  //
  // TwoPRep ITR prefactors
  //

  for (unsigned int p = 0; p < np; p++) {
    char tmp_str[40];
    snprintf(tmp_str, 40, "TwoPRepITR_pfac0_%d", p);
    rdptr vpfac0_ptr(new rdouble(tmp_str));
    TwoPRepITR_vpfac0[p] = vpfac0_ptr;

    snprintf(tmp_str, 40, "TwoPRepITR_pfac1_%d", p);
    rdptr pfac1_ptr(new rdouble(tmp_str));
    TwoPRepITR_pfac1[p] = pfac1_ptr;

    char xyz_str[] = "xyz";
    for (int xyz = 0; xyz < 3; xyz++) {
      char tmp_str[40];
      snprintf(tmp_str, 40, "TwoPRepITR_pfac0_%d_%c", p, xyz_str[xyz]);
      rdptr pfac0_ptr(new rdouble(tmp_str));
      TwoPRepITR_pfac0[p][xyz] = pfac0_ptr;
    }
  }

  //
  // R12_k_G12 VRR prefactors
  //

  rdptr pfac2_ptr(new rdouble("R12kG12_pfac2"));
  R12kG12VRR_pfac2 = pfac2_ptr;

  for (unsigned int p = 0; p < np; p++) {
    char tmp_str[20];
    snprintf(tmp_str, 20, "R12kG12_pfac0_%d", p);
    rdptr vpfac0_ptr(new rdouble(tmp_str));
    R12kG12VRR_vpfac0[p] = vpfac0_ptr;

    snprintf(tmp_str, 20, "R12kG12_pfac1_%d", p);
    rdptr pfac1_ptr(new rdouble(tmp_str));
    R12kG12VRR_pfac1[p] = pfac1_ptr;

    snprintf(tmp_str, 20, "R12kG12_pfac3_%d", p);
    rdptr pfac3_ptr(new rdouble(tmp_str));
    R12kG12VRR_pfac3[p] = pfac3_ptr;

    snprintf(tmp_str, 20, "R12kG12_pfac4_%d", p);
    rdptr vpfac4_ptr(new rdouble(tmp_str));
    R12kG12VRR_vpfac4[p] = vpfac4_ptr;

    char xyz_str[] = "xyz";
    for (int xyz = 0; xyz < 3; xyz++) {
      char tmp_str[20];
      snprintf(tmp_str, 20, "R12kG12_pfac0_%d_%c", p, xyz_str[xyz]);
      rdptr pfac0_ptr(new rdouble(tmp_str));
      R12kG12VRR_pfac0[p][xyz] = pfac0_ptr;

      snprintf(tmp_str, 20, "R12kG12_pfac4_%d_%c", p, xyz_str[xyz]);
      rdptr pfac4_ptr(new rdouble(tmp_str));
      R12kG12VRR_pfac4[p][xyz] = pfac4_ptr;
    }
  }

  for (unsigned int p = 0; p < np; p++) {
    for (int braket = 0; braket < 2; braket++) {
      char tmpstr[200];
      snprintf(tmpstr, 200, "zeta_%c", X[p][braket][0]);
      rdptr zptr(new rdouble(tmpstr));
      zeta[p][braket] = zptr;
    }
  }

  for (unsigned int p = 0; p < np; p++) {
    for (int braket = 0; braket < 2; braket++) {
      char tmpstr[200];
      snprintf(tmpstr, 200, "zeta_%c_2", X[p][braket][0]);
      rdptr zptr(new rdouble(tmpstr));
      zeta2[p][braket] = zptr;
    }
  }

  char alpha12_str[np][20] = {"zeta", "eta"};
  char alpha12_char[np + 1] = "ze";
  for (unsigned int p = 0; p < np; p++) {
    rdptr alpha12_ptr(new rdouble(alpha12_str[p]));
    alpha12[p] = alpha12_ptr;

    char oo2z[20];
    snprintf(oo2z, 20, "oo2%c", alpha12_char[p]);
    rdptr oo2alpha12(new rdouble(oo2z));
    one_o_2alpha12[p] = oo2alpha12;

    char roz[20];
    snprintf(roz, 20, "ro%c", alpha12_char[p]);
    rdptr roz_ptr(new rdouble(roz));
    rho_o_alpha12[p] = roz_ptr;
  }

  //
  // precomputed 1-d ints
  //
  {
    for (unsigned int xyz = 0; xyz != 3; ++xyz) {
      std::ostringstream oss;
      oss << "_0_Overlap_0_" << to_string(CartesianAxis(xyz));
      rdptr p(new rdouble(oss.str()));
      Overlap00_1d[xyz] = p;
    }
  }

#if CTIMEENTITIES_SINGLETONS
  for (unsigned int i = 0; i < NMAX; i++) {
    N_i[i] = prefactor::Scalar((double)i);
  }
#endif
}

Prefactors::~Prefactors() {}

#if CTIMEENTITIES_SINGLETONS
std::shared_ptr<Prefactors::cdouble> Prefactors::Cdouble(double a) {
  return prefactor::Scalar(a);
}
#endif

namespace libint2 {
Prefactors prefactors;
};
