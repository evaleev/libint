/*
 *  Copyright (C) 2004-2024 Edward F. Valeev
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

#include <libint2.hpp>

#include "catch.hpp"

#if defined(LIBINT2_SUPPORT_ERI) && LIBINT2_MAX_AM_eri >= 1

extern "C" {
void init_c_api(unsigned int max_am);

double *compute_eri(unsigned int am1, double alpha1, double *A,
                    unsigned int am2, double alpha2, double *B,
                    unsigned int am3, double alpha3, double *C,
                    unsigned int am4, double alpha4, double *D);

void finalize_c_api();
}

TEST_CASE("C API", "[c-api]") {
  int am1, am2, am3, am4;
  double alpha1, alpha2, alpha3, alpha4;
  double A[3], B[3], C[3], D[3];
  am1 = am2 = am3 = am4 = 1;
  alpha1 = 1.1;
  alpha2 = 2.3;
  alpha3 = 3.4;
  alpha4 = 4.8;
  A[0] = 0.0;
  A[1] = 1.0;
  A[2] = 2.0;
  B[0] = 1.0;
  B[1] = 2.0;
  B[2] = 0.0;
  C[0] = 2.0;
  C[1] = 0.0;
  C[2] = 1.0;
  D[0] = 0.0;
  D[1] = 1.0;
  D[2] = 2.0;

  using std::max;
  auto max_am = max(max(am1, am2), max(am3, am4));
  init_c_api(max_am);

  auto *c_result = compute_eri(am1, alpha1, A, am2, alpha2, B, am3, alpha3, C,
                               am4, alpha4, D);

  const double *cpp_result;
  using libint2::Shell;
  Shell sh1{{alpha1}, {{am1, false, {1.0}}}, {A[0], A[1], A[2]}};
  Shell sh2{{alpha2}, {{am2, false, {1.0}}}, {B[0], B[1], B[2]}};
  Shell sh3{{alpha3}, {{am3, false, {1.0}}}, {C[0], C[1], C[2]}};
  Shell sh4{{alpha4}, {{am4, false, {1.0}}}, {D[0], D[1], D[2]}};
  libint2::Engine engine(libint2::Operator::coulomb, 1, max_am);
  engine.compute(sh1, sh2, sh3, sh4);
  cpp_result = engine.results()[0];

  printf("CMake Configuration (C)  : %s\n", configuration_accessor());
  printf("CMake Configuration (C++): %s\n",
         libint2::configuration_accessor().c_str());

  unsigned int n1, n2, n3, n4;
  int a, b, c, d, abcd;
  n1 = (am1 + 1) * (am1 + 2) / 2;
  n2 = (am2 + 1) * (am2 + 2) / 2;
  n3 = (am3 + 1) * (am3 + 2) / 2;
  n4 = (am4 + 1) * (am4 + 2) / 2;
  const auto norm_factor = sh1.contr[0].coeff[0] * sh2.contr[0].coeff[0] *
                           sh3.contr[0].coeff[0] * sh4.contr[0].coeff[0];
  for (a = 0, abcd = 0; a < n1; a++) {
    for (b = 0; b < n2; b++) {
      for (c = 0; c < n3; c++) {
        for (d = 0; d < n4; d++, ++abcd) {
          // printf("a = %d b = %d c = %d d = %d (ab|cd) = %20.15lf , ref
          // (ab|cd) = %20.15lf\n", a, b, c, d, c_result[abcd]*norm_factor,
          // cpp_result[abcd]);
          REQUIRE(c_result[abcd] * norm_factor == Approx(cpp_result[abcd]));
        }
      }
    }
  }

  finalize_c_api();
}

#endif
