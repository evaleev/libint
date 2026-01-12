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

#include "catch.hpp"
#include "fixture.h"

#define LIBINT2_TEST_ONEBODY(scale, op_idx)                         \
  int m0, m1;                                                       \
  FOR_SOLIDHARM_STANDARD(l0, m0)                                    \
  FOR_SOLIDHARM_STANDARD(l1, m1)                                    \
  const auto i0_s = libint2::INT_SOLIDHARMINDEX_STANDARD(l0, m0);   \
  const auto i1_s = libint2::INT_SOLIDHARMINDEX_STANDARD(l1, m1);   \
  const auto i0i1_s = i0_s * n1 + i1_s;                             \
  if (libint2::solid_harmonics_ordering() ==                        \
      libint2::SHGShellOrdering_Standard) {                         \
    REQUIRE(engine.results()[op_idx][i0i1_s] / scale ==             \
            Approx(shellset_ref_standard[i0i1_s]));                 \
  } else {                                                          \
    const auto i0_g = libint2::INT_SOLIDHARMINDEX_GAUSSIAN(l0, m0); \
    const auto i1_g = libint2::INT_SOLIDHARMINDEX_GAUSSIAN(l1, m1); \
    const auto i0i1_g = i0_g * n1 + i1_g;                           \
    REQUIRE(engine.results()[op_idx][i0i1_g] / scale ==             \
            Approx(shellset_ref_standard[i0i1_s]));                 \
  }                                                                 \
  END_FOR_SOLIDHARM                                                 \
  END_FOR_SOLIDHARM

TEST_CASE_METHOD(libint2::unit::DefaultFixture, "electrostatic potential",
                 "[engine][1-body]") {
#if defined(LIBINT2_SUPPORT_ONEBODY)

  std::vector<Shell> obs{
      Shell{{1.0, 3.0}, {{2, true, {1.0, 0.3}}}, {{0.0, 0.0, 0.0}}},
      Shell{{2.0, 5.0}, {{2, true, {1.0, 0.2}}}, {{1.0, 1.0, 1.0}}}};

  {
    const auto lmax = std::min(3, LIBINT2_MAX_AM_elecpot);
    if (lmax >= 2) {
      auto engine = Engine(Operator::nuclear, 2, lmax);
      engine.set_params(make_point_charges(atoms));

      const auto scale = 2.3;
      engine.prescale_by(scale);
      engine.compute(obs[0], obs[0]);
      {
        const auto l0 = obs[0].contr[0].l;
        const auto l1 = obs[0].contr[0].l;
        const auto n1 = 2 * l1 + 1;

        // this is laid out in standard solids order
        std::vector<double> shellset_ref_standard = {
            -1.238239259091998e+01, 0.000000000000000e+00,
            0.000000000000000e+00,  -5.775996163160049e-02,
            0.000000000000000e+00,  0.000000000000000e+00,
            -1.301230978657952e+01, -6.796143730068988e-02,
            0.000000000000000e+00,  1.139389632827834e-01,
            0.000000000000000e+00,  -6.796143730068988e-02,
            -1.343732979153083e+01, 0.000000000000000e+00,
            -1.478824785355970e-02, -5.775996163160049e-02,
            0.000000000000000e+00,  0.000000000000000e+00,
            -1.284475452992947e+01, 0.000000000000000e+00,
            0.000000000000000e+00,  1.139389632827834e-01,
            -1.478824785355970e-02, 0.000000000000000e+00,
            -1.241040347301479e+01};

        LIBINT2_TEST_ONEBODY(scale, 0);
      }

      engine.prescale_by(1);
      engine.compute(obs[0], obs[1]);
      {
        const auto l0 = obs[0].contr[0].l;
        const auto l1 = obs[1].contr[0].l;
        const auto n1 = 2 * l1 + 1;

        // this is laid out in standard solids order
        std::vector<double> shellset_ref_standard = {
            -4.769186621041819e-01, -9.303619356400431e-01,
            -1.559058302243514e+00, -9.290824121864600e-01,
            -5.835786921473129e-04, -1.159266418436018e+00,
            -3.770080831197964e-01, 9.572841308198474e-01,
            -8.291498398421207e-01, -1.663667687168316e+00,
            -2.171951144148577e+00, 1.074249956874296e+00,
            2.128355904665372e+00,  1.074590109905394e+00,
            -3.485163651594458e-03, -1.160865205880651e+00,
            -8.344173649626901e-01, 9.566621490332916e-01,
            -3.760919234260182e-01, 1.660514988916377e+00,
            -1.120272634615116e-03, -1.385603731947886e+00,
            -2.105750177166632e-03, 1.380654897976564e+00,
            2.115041199099945e+00};

        LIBINT2_TEST_ONEBODY(1.0, 0);
      }
    }
  }

  // see https://github.com/evaleev/libint/issues/199
  {
    const auto lmax = std::min(3, LIBINT2_MAX_AM_elecpot);
    if (lmax >= 2) {
      const auto deriv_order = INCLUDE_ONEBODY;
      auto engine = Engine(Operator::nuclear, 2, lmax, deriv_order);
      engine.set_params(make_point_charges(atoms));
      const auto& buf = engine.results();
      REQUIRE(libint2::num_geometrical_derivatives(atoms.size() + 2,
                                                   deriv_order) == buf.size());
    }
  }

#endif  // LIBINT2_SUPPORT_ONEBODY
}

TEST_CASE_METHOD(libint2::unit::DefaultFixture, "erf correctness",
                 "[engine][1-body]") {
#if defined(LIBINT2_SUPPORT_ONEBODY)
  if (LIBINT_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD) return;

  constexpr int l0 = 2;
  constexpr int l1 = 2;
  constexpr int n1 = 2 * l1 + 1;
  std::vector<Shell> obs{
      Shell{{1.0, 3.0}, {{l0, true, {1.0, 0.3}}}, {{0.0, 0.0, 0.0}}},
      Shell{{2.0, 5.0}, {{l1, true, {1.0, 0.2}}}, {{1.0, 1.0, 1.0}}}};
  {
    const auto lmax = std::min(3, LIBINT2_MAX_AM_elecpot);
    if (lmax >= 2) {
      auto engine = Engine(Operator::nuclear, 2, lmax);
      engine.set_params(make_point_charges(atoms));

      engine.compute(obs[0], obs[1]);
      {
        // this is laid out in standard solids order
        std::vector<double> shellset_ref_standard = {
            -4.769186621041819e-01, -9.303619356400431e-01,
            -1.559058302243514e+00, -9.290824121864600e-01,
            -5.835786921473129e-04, -1.159266418436018e+00,
            -3.770080831197964e-01, 9.572841308198474e-01,
            -8.291498398421207e-01, -1.663667687168316e+00,
            -2.171951144148577e+00, 1.074249956874296e+00,
            2.128355904665372e+00,  1.074590109905394e+00,
            -3.485163651594458e-03, -1.160865205880651e+00,
            -8.344173649626901e-01, 9.566621490332916e-01,
            -3.760919234260182e-01, 1.660514988916377e+00,
            -1.120272634615116e-03, -1.385603731947886e+00,
            -2.105750177166632e-03, 1.380654897976564e+00,
            2.115041199099945e+00};
        LIBINT2_TEST_ONEBODY(1.0, 0);
      }
    }
  }
#endif  // LIBINT2_SUPPORT_ONEBODY
}

TEST_CASE_METHOD(libint2::unit::DefaultFixture, "W correctness",
                 "[engine][1-body]") {
#if defined(LIBINT2_SUPPORT_ONEBODY)
  if (LIBINT_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD) return;

  constexpr int l0 = 2;
  constexpr int l1 = 2;
  constexpr int n1 = 2 * l1 + 1;
  std::vector<Shell> obs{
      Shell{{1.0, 3.0}, {{l0, true, {1.0, 0.3}}}, {{0.0, 0.0, 0.0}}},
      Shell{{2.0, 5.0}, {{l1, true, {1.0, 0.2}}}, {{1.0, 1.0, 1.0}}}};
  {
    const auto lmax = std::min(3, LIBINT2_MAX_AM_elecpot);
    if (lmax >= 2) {
      auto engine = Engine(Operator::opVop, 2, lmax);
      engine.set_params(make_point_charges(atoms));

      engine.compute(obs[0], obs[1]);
      // all ref are laid out in standard solids order
      {
        std::vector<double> shellset_ref_standard = {
            -8.01340642466483e+00, -1.30377852807092e+01, -9.93488915362642e+00,
            -1.30205110577182e+01, 1.49417100784088e-03,  -1.53323266823052e+01,
            -6.93468224743044e+00, 7.09818366817654e+00,  -1.26718835918618e+01,
            -1.05956288853455e+01, -1.32056241932669e+01, 5.91008766684120e+00,
            1.64117141263422e+01,  5.91170569953823e+00,  -9.97707886886801e-03,
            -1.53405810851944e+01, -1.27006169503565e+01, 7.08618624960792e+00,
            -6.92922348322532e+00, 1.05765790465850e+01,  -1.01288023970354e-02,
            -9.37479473854916e+00, -8.12682563116951e-03, 9.35287687372449e+00,
            1.66715357961274e+01};
        LIBINT2_TEST_ONEBODY(1.0, 0);
      }
      {
        std::vector<double> shellset_ref_standard = {
            -7.88373054495278e-01, -5.92741140473083e-01, 7.18456828001176e-02,
            -6.39138844263037e-01, 2.20726461230564e-01,  -5.92585516679938e-01,
            2.34732079397513e-01,  1.73416688435562e+00,  2.02375318329227e+00,
            4.29028634252541e-02,  6.08798471053324e-01,  1.42709886052571e+00,
            -3.62051356103073e-01, 2.58807312093735e+00,  6.88208626980571e-01,
            1.62736755338176e-01,  1.16308344349130e+00,  9.78124133242033e-01,
            1.15541165536756e+00,  -1.13849207132631e-02, 4.18505962024446e-01,
            1.90312514813159e-01,  -6.40076159142129e-01, 8.81020773424062e-01,
            1.05645673339062e+00};
        LIBINT2_TEST_ONEBODY(1.0, 1);
      }
      {
        std::vector<double> shellset_ref_standard = {
            7.89010735420826e-01,  6.29303049284450e-01,  -6.43094112261234e-02,
            5.72123450669670e-01,  1.98231144425269e-01,  -1.52481061826136e-01,
            -1.15201706049790e+00, -9.59578826810163e-01, -1.15664268784483e+00,
            2.53464418195577e-04,  -6.01303780133668e-01, -2.58645514834283e+00,
            3.75234373832125e-01,  -1.43968950861751e+00, 7.18031641380203e-01,
            6.13802108561365e-01,  -2.03452230861181e+00, -1.70486591152126e+00,
            -2.38135136433799e-01, 3.34943994067535e-02,  4.38398996368850e-01,
            8.80143879113653e-01,  -6.47475279081898e-01, 1.99006135765209e-01,
            -1.06105366740987e+00};
        LIBINT2_TEST_ONEBODY(1.0, 2);
      }
      {
        std::vector<double> shellset_ref_standard = {
            -2.49761988715857e-04, 1.65958345836693e+00,  3.34242295115910e-03,
            -1.66148554416249e+00, -3.03391179873709e+00, 1.49079297243232e+00,
            8.00331729579887e-01,  -6.80731543814720e-01, 4.72634239043335e-01,
            -1.19318769283081e+00, 1.11521084365706e-03,  1.33833228953250e-02,
            -6.62037517969870e-03, 1.19969136764911e-02,  -1.05965708465147e-03,
            -1.48836792336934e+00, -4.67132013887958e-01, 6.74189183691556e-01,
            -8.03426805711286e-01, -1.18966404475986e+00, -1.44170566149487e+00,
            -2.10343679084144e+00, -1.89380199837230e+00, -2.10786126986013e+00,
            7.49910831098100e-04};
        LIBINT2_TEST_ONEBODY(1.0, 3);
      }
    }
  }
#endif  // LIBINT2_SUPPORT_ONEBODY
}

// verify that python/tests/test_libint2.py:test_integrals is correct
TEST_CASE_METHOD(libint2::unit::DefaultFixture, "python correctness",
                 "[engine][1-body]") {
#if defined(LIBINT2_SUPPORT_ONEBODY)
  std::vector<Shell> obs{Shell{{1.0}, {{1, true, {10.0}}}, {{0.0, 0.0, 0.0}}},
                         Shell{{1.0}, {{2, true, {10.0}}}, {{0.1, 0.2, 0.3}}}};
  constexpr int l0 = 1;
  constexpr int l1 = 2;
  constexpr int n1 = 2 * l1 + 1;

  {
    const auto lmax = LIBINT2_MAX_AM_overlap;
    if (lmax >= 2) {
      auto engine = Engine(Operator::overlap, 2, lmax);

      engine.compute(obs[0], obs[1]);

      // this is laid out in standard solids order
      //      REQUIRE(engine.results()[0][0] == Approx(-0.08950980671097111));
      //      REQUIRE(engine.results()[0][1] == Approx(-0.2685294201329133));
      //      REQUIRE(engine.results()[0][5] == Approx(0.0055943629194356937));
      std::vector<double> shellset_ref_standard = {
          -0.0895098067109711,  -0.268529420132913, 0.114661696280563,
          0.00559436291943569,  0.183681582521472,  0.00559436291943569,
          -0.169695675222883,   -0.312493496201254, -0.0848478376114413,
          -0.00419577218957676, -0.184613976341378, 0.00559436291943569,
          0.0573308481402817,   -0.276920964512067, -0.0946379727204538};

      LIBINT2_TEST_ONEBODY(1.0, 0);
    }
  }
#endif  // LIBINT2_SUPPORT_ONEBODY
}
