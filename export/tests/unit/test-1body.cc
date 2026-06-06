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
      const auto deriv_order = LIBINT_INCLUDE_ONEBODY;
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

TEST_CASE_METHOD(libint2::unit::DefaultFixture, "σpeμσp correctness",
                 "[engine][1-body]") {
#if defined(LIBINT2_SUPPORT_ONEBODY)
  if (LIBINT_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD) return;

  // Two contracted Gaussian shells at distinct centers. We exercise
  // `(σ·p) r (σ·p)` over the (s|d) and (d|s) shell pairs and validate
  // Hermiticity: trace components (q=0) are symmetric under bra↔ket swap,
  // antisym components (q=1,2,3) are antisymmetric.
  std::vector<Shell> obs{
      Shell{{1.0, 3.0}, {{0, true, {1.0, 0.3}}}, {{0.0, 0.0, 0.0}}},
      Shell{{2.0, 5.0}, {{2, true, {1.0, 0.2}}}, {{1.0, 1.0, 1.0}}}};

  const auto lmax = std::min(2, LIBINT2_MAX_AM_opemuop);
  if (lmax < 2) return;

  auto engine = Engine(Operator::σpeμσp, 2, lmax);
  engine.set_params(std::array<double, 3>{{0.0, 0.0, 0.0}});

  // (s|σpeμσp|d) and (d|σpeμσp|s)
  engine.compute(obs[0], obs[1]);
  std::array<std::vector<double>, 12> ab;
  for (int c = 0; c < 12; ++c) {
    const auto* buf = engine.results()[c];
    REQUIRE(buf != nullptr);
    ab[c].assign(buf, buf + (1 * 5));  // n_s × n_d_pure = 1 × 5
  }

  engine.compute(obs[1], obs[0]);
  std::array<std::vector<double>, 12> ba;
  for (int c = 0; c < 12; ++c) {
    const auto* buf = engine.results()[c];
    REQUIRE(buf != nullptr);
    ba[c].assign(buf, buf + (5 * 1));  // n_d_pure × n_s
  }

  // Hermiticity check: σpeμσp is Hermitian, but the Pauli identity routes the
  // imaginary i factor on the antisym pieces into a real-stored sign flip.
  //   q=0 trace: matrix is symmetric  ⇒  ab[0+k][i,j] ==  ba[0+k][j,i]
  //   q=1..3   : matrix is antisym    ⇒  ab[q+k][i,j] == -ba[q+k][j,i]
  // (Indices: ab is laid out (i=0..n_s-1, j=0..n_d-1); ba is (j, i).)
  const double tol = 1.0e-12;
  for (int k = 0; k < 3; ++k) {
    for (int q = 0; q < 4; ++q) {
      const int comp = 4 * k + q;
      const double expected_sign = (q == 0) ? +1.0 : -1.0;
      for (int i = 0; i < 1; ++i) {
        for (int j = 0; j < 5; ++j) {
          const double v_ab = ab[comp][i * 5 + j];
          const double v_ba = ba[comp][j * 1 + i];
          REQUIRE(std::isfinite(v_ab));
          REQUIRE(std::abs(v_ab - expected_sign * v_ba) < tol);
        }
      }
    }
  }

  // Sanity: not every component is identically zero (would mask codegen bugs).
  bool any_nonzero = false;
  for (int c = 0; c < 12; ++c)
    for (double v : ab[c])
      if (std::abs(v) > 1.0e-10) any_nonzero = true;
  REQUIRE(any_nonzero);
#endif  // LIBINT2_SUPPORT_ONEBODY
}


// Helper: max number of primitives across all centers in q_gau data.
static size_t max_nprim(const libint2::GaussianPotentialCentersData& data) {
  size_t result = 0;
  for (const auto& ptr : data)
    if (ptr) result = std::max(result, ptr->size());
  return result;
}

// Helper: compute lower-triangle (packed) 1-body matrix for a nuclear-type
// operator using the given engine. Returns n*(n+1)/2 elements.
static std::vector<double> compute_nuclear_ltri(Engine& engine,
                                                const BasisSet& obs) {
  const auto n = libint2::nbf(obs);
  const auto nshells = obs.size();
  auto shell2bf = obs.shell2bf();
  const auto& buf = engine.results();
  const auto ntri = n * (n + 1) / 2;
  std::vector<double> V(ntri, 0.0);
  for (size_t s1 = 0; s1 < nshells; ++s1) {
    auto bf1 = shell2bf[s1];
    auto n1 = obs[s1].size();
    for (size_t s2 = 0; s2 <= s1; ++s2) {
      auto bf2 = shell2bf[s2];
      auto n2 = obs[s2].size();
      engine.compute(obs[s1], obs[s2]);
      for (size_t i = 0; i < n1; ++i)
        for (size_t j = 0; j < n2; ++j) {
          const auto ii = bf1 + i;
          const auto jj = bf2 + j;
          if (ii >= jj) V[ii * (ii + 1) / 2 + jj] = buf[0][i * n2 + j];
        }
    }
  }
  return V;
}

// Helper: compare two lower-triangle vectors element-wise.
static void compare_ltri(const std::vector<double>& V,
                         const std::vector<double>& V_ref, size_t n, double eps,
                         const std::string& label) {
  REQUIRE(V.size() == V_ref.size());
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      const auto idx = i * (i + 1) / 2 + j;
      INFO(label << "(" << i << "," << j << ") = " << V[idx]
                 << " ref = " << V_ref[idx]);
      if (std::abs(V_ref[idx]) > 1e-12) {
        REQUIRE(V[idx] == Approx(V_ref[idx]).epsilon(eps));
      } else {
        REQUIRE(std::abs(V[idx]) < std::max(eps, 1e-14));
      }
    }
  }
}

TEST_CASE_METHOD(libint2::unit::DefaultFixture,
                 "q_gau point nuclear matches Operator::nuclear",
                 "[engine][1-body]") {
#if defined(LIBINT2_SUPPORT_ONEBODY)
  BasisSet obs("sto-3g", atoms);
  auto point_charges = make_point_charges(atoms);
  const auto n = libint2::nbf(obs);

  // Reference: Operator::nuclear
  Engine nuc_engine(Operator::nuclear, obs.max_nprim(), obs.max_l());
  nuc_engine.set_params(point_charges);
  auto V_nuc = compute_nuclear_ltri(nuc_engine, obs);

  // q_gau with point nuclear model
  auto q_gau_data =
      libint2::make_q_gau_data(libint2::NuclearModel::PointCharge, atoms);
  Engine q_engine(Operator::q_gau,
                  std::max(obs.max_nprim(), max_nprim(q_gau_data)),
                  obs.max_l());
  q_engine.set_params(std::make_tuple(q_gau_data, point_charges));
  auto V_qgau = compute_nuclear_ltri(q_engine, obs);

  compare_ltri(V_qgau, V_nuc, n, 1e-14, "q_gau_pt");
#endif
}

TEST_CASE_METHOD(libint2::unit::DefaultFixture,
                 "q_gau erf matches Operator::erf_nuclear",
                 "[engine][1-body]") {
#if defined(LIBINT2_SUPPORT_ONEBODY)
  BasisSet obs("sto-3g", atoms);
  auto point_charges = make_point_charges(atoms);
  const auto n = libint2::nbf(obs);
  const double omega = 0.5;

  // Reference: Operator::erf_nuclear
  Engine erf_engine(Operator::erf_nuclear, obs.max_nprim(), obs.max_l());
  erf_engine.set_params(std::make_tuple(omega, point_charges));
  auto V_erf = compute_nuclear_ltri(erf_engine, obs);

  // q_gau with erf model
  auto q_gau_data = libint2::make_q_gau_data_erf(omega, atoms);
  Engine q_engine(Operator::q_gau,
                  std::max(obs.max_nprim(), max_nprim(q_gau_data)),
                  obs.max_l());
  q_engine.set_params(std::make_tuple(q_gau_data, point_charges));
  auto V_qgau = compute_nuclear_ltri(q_engine, obs);

  compare_ltri(V_qgau, V_erf, n, 1e-14, "q_gau_erf");
#endif
}

TEST_CASE_METHOD(libint2::unit::DefaultFixture,
                 "q_gau erfc matches Operator::erfc_nuclear",
                 "[engine][1-body]") {
#if defined(LIBINT2_SUPPORT_ONEBODY)
  BasisSet obs("sto-3g", atoms);
  auto point_charges = make_point_charges(atoms);
  const auto n = libint2::nbf(obs);
  const double omega = 0.5;

  // Reference: Operator::erfc_nuclear
  Engine erfc_engine(Operator::erfc_nuclear, obs.max_nprim(), obs.max_l());
  erfc_engine.set_params(std::make_tuple(omega, point_charges));
  auto V_erfc = compute_nuclear_ltri(erfc_engine, obs);

  // q_gau with erfc model
  auto q_gau_data = libint2::make_q_gau_data_erfc(omega, atoms);
  Engine q_engine(Operator::q_gau,
                  std::max(obs.max_nprim(), max_nprim(q_gau_data)),
                  obs.max_l());
  q_engine.set_params(std::make_tuple(q_gau_data, point_charges));
  auto V_qgau = compute_nuclear_ltri(q_engine, obs);

  compare_ltri(V_qgau, V_erfc, n, 1e-14, "q_gau_erfc");
#endif
}

TEST_CASE_METHOD(libint2::unit::DefaultFixture,
                 "q_gau erfx matches Operator::erfx_nuclear",
                 "[engine][1-body]") {
#if defined(LIBINT2_SUPPORT_ONEBODY)
  BasisSet obs("sto-3g", atoms);
  auto point_charges = make_point_charges(atoms);
  const auto n = libint2::nbf(obs);
  const double omega = 0.5, lambda = 0.3, sigma = 0.7;

  // Reference: Operator::erfx_nuclear
  Engine erfx_engine(Operator::erfx_nuclear, obs.max_nprim(), obs.max_l());
  erfx_engine.set_params(std::make_tuple(
      std::array<double, 3>{omega, lambda, sigma}, point_charges));
  auto V_erfx = compute_nuclear_ltri(erfx_engine, obs);

  // q_gau with erfx model
  auto q_gau_data = libint2::make_q_gau_data_erfx(omega, lambda, sigma, atoms);
  Engine q_engine(Operator::q_gau,
                  std::max(obs.max_nprim(), max_nprim(q_gau_data)),
                  obs.max_l());
  q_engine.set_params(std::make_tuple(q_gau_data, point_charges));
  auto V_qgau = compute_nuclear_ltri(q_engine, obs);

  compare_ltri(V_qgau, V_erfx, n, 1e-14, "q_gau_erfx");
#endif
}

TEST_CASE_METHOD(libint2::unit::DefaultFixture,
                 "q_gau Gaussian nuclear matches per-center erf_nuclear",
                 "[engine][1-body]") {
#if defined(LIBINT2_SUPPORT_ONEBODY)
  // Gaussian nuclear model: V(r) = -(eZ/r)*erf(sqrt(xi)*r)
  // equivalent to erf_nuclear with omega = sqrt(xi(Z)) per center.
  // For a heteronuclear system, build reference by accumulating per-center
  // erf_nuclear contributions, each with omega = sqrt(xi(Z_i)).
  BasisSet obs("sto-3g", atoms);
  auto point_charges = make_point_charges(atoms);
  const auto n = libint2::nbf(obs);
  const auto ntri = n * (n + 1) / 2;

  // Reference: accumulate per-center erf_nuclear with omega = sqrt(xi(Z))
  std::vector<double> V_ref(ntri, 0.0);
  for (size_t a = 0; a < atoms.size(); ++a) {
    const int Z = atoms[a].atomic_number;
    if (Z == 0) continue;
    const double xi = libint2::chemistry::gaussian_nuclear_exponent(Z);
    const double omega = std::sqrt(xi);
    // single-center point charge
    std::vector<std::pair<double, std::array<double, 3>>> single_charge = {
        point_charges[a]};
    Engine erf_engine(Operator::erf_nuclear, obs.max_nprim(), obs.max_l());
    erf_engine.set_params(std::make_tuple(omega, single_charge));
    auto V_center = compute_nuclear_ltri(erf_engine, obs);
    for (size_t k = 0; k < ntri; ++k) V_ref[k] += V_center[k];
  }

  // q_gau with Gaussian nuclear model
  auto q_gau_data =
      libint2::make_q_gau_data(libint2::NuclearModel::GaussianCharge, atoms);
  Engine q_engine(Operator::q_gau,
                  std::max(obs.max_nprim(), max_nprim(q_gau_data)),
                  obs.max_l());
  q_engine.set_params(std::make_tuple(q_gau_data, point_charges));
  auto V_qgau = compute_nuclear_ltri(q_engine, obs);

  compare_ltri(V_qgau, V_ref, n, 1e-14, "q_gau_finite_nuc");
#endif
}

TEST_CASE_METHOD(libint2::unit::DefaultFixture, "q_gau SAP correctness",
                 "[engine][1-body]") {
#if defined(LIBINT2_SUPPORT_ONEBODY)
  // atoms from fixture: O(0,0,0), O(0,0,2), H(0,-1,-1), H(0,1,3) in Bohr
  BasisSet obs("sto-3g", atoms);
  auto point_charges = make_point_charges(atoms);
  const auto n = libint2::nbf(obs);

  auto q_gau_data = libint2::make_q_gau_data(libint2::NuclearModel::PointCharge,
                                             atoms, "sap_helfem_large");
  Engine q_engine(Operator::q_gau,
                  std::max(obs.max_nprim(), max_nprim(q_gau_data)),
                  obs.max_l());
  q_engine.set_params(std::make_tuple(q_gau_data, point_charges));
  auto V_sap = compute_nuclear_ltri(q_engine, obs);

  // clang-format off
  // Reference SAP lower triangle (packed, 78 elements for 12x12 matrix).
  // See source directory for reference generation details.
  const double V_sap_ref[] = {
      -48.161626712664187,
      -4.7836856246433701, -3.5844048783736095,
      0, 0, -3.3343483508728493,
      0.012917426701848284, 0.13411995427355725, 0, -3.4287357703538213,
      -0.010967831892154056, -0.23517850034936605, 0, -0.096970266316132495,
      -3.7771303420092082,
      -2.0033821099010749e-06, -0.77175616234890021, 0,
      -0.00034296744926888052, -1.2810294579879118, -48.161626712664173,
      -0.77175616234890032, -1.2542321767214553, 0, 0.0063225004608634238,
      -1.378150375719013, -4.7836856246433701, -3.5844048783736095,
      0, 0, -0.48800756489751196, 0, 0, 0, 0, -3.3343483508728493,
      0.00034296744926888063, -0.0063225004608634264, 0,
      -0.49670368376476004, -0.02113524909099768, -0.012917426701848284,
      -0.13411995427355725, 0, -3.4287357703538213,
      1.2810294579879122, 1.3781503757190121, 0, -0.021135249090997673,
      1.1856000889763396, 0.010967831892154431, 0.23517850034936605, 0,
      -0.096970266316132495, -3.7771303420092082,
      -1.8612853433654166, -1.7121869000279206, 0, 0.79173492433793191,
      0.66883426856379913, -0.21984583270685026, -0.40830320535018161, 0,
      0.076419612610196419, 0.43074554813761501, -1.8237300358596116,
      -0.21984583270685035, -0.40830320535018183, 0, -0.07641961261019653,
      -0.43074554813761523, -1.8612853433654166, -1.7121869000279188, 0,
      -0.79173492433793236, -0.66883426856379957, -0.14206535846913848,
      -1.8237300358596116};
  // clang-format on

  const auto ntri = n * (n + 1) / 2;
  std::vector<double> V_ref(V_sap_ref, V_sap_ref + ntri);
  compare_ltri(V_sap, V_ref, n, 1e-10, "q_gau_sap");
#endif
}

// Helper: compute full n×n 1-body matrix for a specific component of a
// multi-component nuclear-type operator, then pack as lower triangle.
static std::vector<double> compute_nuclear_ltri_comp(Engine& engine,
                                                     const BasisSet& obs,
                                                     int comp) {
  const auto n = libint2::nbf(obs);
  const auto nshells = obs.size();
  auto shell2bf = obs.shell2bf();
  const auto& buf = engine.results();
  const auto ntri = n * (n + 1) / 2;
  // Compute full matrix (need all shell pairs for accumulation over centers)
  std::vector<double> M(n * n, 0.0);
  for (size_t s1 = 0; s1 < nshells; ++s1) {
    auto bf1 = shell2bf[s1];
    auto n1 = obs[s1].size();
    for (size_t s2 = 0; s2 < nshells; ++s2) {
      auto bf2 = shell2bf[s2];
      auto n2 = obs[s2].size();
      engine.compute(obs[s1], obs[s2]);
      if (buf[comp] == nullptr) continue;
      for (size_t i = 0; i < n1; ++i)
        for (size_t j = 0; j < n2; ++j)
          M[(bf1 + i) * n + (bf2 + j)] = buf[comp][i * n2 + j];
    }
  }
  // Pack lower triangle
  std::vector<double> V(ntri, 0.0);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j <= i; ++j) V[i * (i + 1) / 2 + j] = M[i * n + j];
  return V;
}

TEST_CASE_METHOD(libint2::unit::DefaultFixture, "op_q_gau_op SAP correctness",
                 "[engine][1-body]") {
#if defined(LIBINT2_SUPPORT_ONEBODY)
  // atoms from fixture: O(0,0,0), O(0,0,2), H(0,-1,-1), H(0,1,3) in Bohr
  // sto-3g gives n=12 basis functions
  BasisSet obs("sto-3g", atoms);
  auto point_charges = make_point_charges(atoms);
  const auto n = libint2::nbf(obs);

  // Build SAP-only data: replace the point nuclear {inf, 1.0} with {inf, 0.0}
  // so only the SAP correction primitives contribute (no bare Coulomb).
  auto q_gau_data = libint2::make_q_gau_data(libint2::NuclearModel::PointCharge,
                                             atoms, "sap_grasp_large");
  libint2::GaussianPotentialCentersData sap_only_data;
  sap_only_data.reserve(q_gau_data.size());
  for (const auto& ptr : q_gau_data) {
    if (ptr && !ptr->empty()) {
      auto data = *ptr;           // copy
      data[0].coefficient = 0.0;  // zero out the point nuclear coefficient
      sap_only_data.push_back(
          std::make_shared<const libint2::GaussianPotentialData>(
              std::move(data)));
    } else {
      sap_only_data.push_back(ptr);
    }
  }
  const auto lmax =
      std::min(static_cast<int>(obs.max_l()), LIBINT2_MAX_AM_elecpot);
  Engine q_engine(Operator::op_q_gau_op,
                  std::max(obs.max_nprim(), max_nprim(sap_only_data)), lmax);
  q_engine.set_params(std::make_tuple(sap_only_data, point_charges));

  // clang-format off
  // Reference W_sap lower triangles from external 4c Dirac code (real parts).
  // These are the SAP correction only (without the bare Coulomb opVop).
  // Component 0 (scalar), 78 elements
  const double W_ref_s[] = {
      1017.1940316588756,
      -6.7517438914914019, 14.90263333586061,
      0, 0, 59.231740315770757,
      0.0067660273198012749, -0.032326959926133783, 0, 59.238689218930077,
      -0.22936483274329972, 1.6831735948256663, 0, 0.0096263636335646371, 59.881035733738578,
      -0.00094922461392558217, -0.49401374578031587, 0, -0.0010286720297784421, -0.74258885324890422, 1017.1940316588756,
      -0.49401374577993529, 1.0557146456752298, 0, -0.0019906557365285044, 5.3008728671145358, -6.7517438914914036, 14.902633335860607,
      0, 0, 2.972026321615155, 0, 0, 0, 0, 59.231740315770757,
      0.0010286720297784711, 0.0019906557365284975, 0, 2.9711640731200957, 0.001938293450133753, -0.0067660273198012887, 0.032326959926133783, 0, 59.238689218930077,
      0.74258885324854762, -5.3008728671145224, 0, 0.0019382934501337695, -12.733566048947557, 0.22936483274331063, -1.6831735948256625, 0, 0.0096263636335646476, 59.881035733738578,
      0.68020575472394618, 3.8324410553887276, 0, -6.1334779485505644, -5.7614548706357862, -0.098493544021966237, -0.5349072985870249, 0, -0.096640899779031425, -0.025753311617736457, 11.034319903150422,
      -0.098493544021899443, -0.53490729858702546, 0, 0.096640899779031925, 0.025753311617739822, 0.68020575472423772, 3.8324410553887276, 0, 6.1334779485505653, 5.7614548706357844, -0.30378068899204275, 11.034319903150422};
  // Component 1 (x), 78 elements
  const double W_ref_x[] = {
      0,
      0, 0,
      0, 0, 0,
      -1.5576470863370833, -0.90699161564122788, 0, 0,
      -0.04151583448537164, -0.018843269537484188, 0, -6.7719625920746838, 0,
      0, 0.0022578738355661931, 0, -0.057002833209069528, 0.0036697623210374501, 0,
      -0.0022578738355661949, 0, 0, -0.97547564303023371, 0.0050498697820343993, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      -0.057002833209069965, -0.9754756430302336, 0, 0, -0.83353261538165768, 1.5576470863370886, 0.9069916156412271, 0, 0,
      0.0036697623210374588, 0.0050498697820344045, 0, 0.83353261538165646, 0, 0.041515834485371626, 0.018843269537484185, 0, -6.7719625920746838, 0,
      0.10694989397838189, 0.21179524040308861, 0, 1.2532746958054455, -1.1391504297329242, -0.0071399526451828427, -0.14271807794149374, 0, 0.34084650640231423, 0.094452927305536089, 0,
      -0.0071399526451828427, -0.14271807794149355, 0, -0.34084650640231434, -0.094452927305535311, 0.10694989397838202, 0.21179524040308961, 0, -1.2532746958054461, 1.1391504297329258, 0, 0};
  // Component 2 (y), 78 elements
  const double W_ref_y[] = {
      0,
      0, 0,
      1.5576470863370833, 0.90699161564122788, 0,
      0, 0, -0.013127348731989879, 0,
      0, 0, 6.7819650765550188, 0, 0,
      0, 0, 0.057099410762780578, 0, 0, 0,
      0, 0, 0.97887720985746873, 0, 0, 0, 0,
      0.057099410762781022, 0.97887720985746862, 0, 0.0033118181997882313, 0.83646315449674191, -1.5576470863370886, -0.9069916156412271, 0,
      0, 0, -0.0033118181997882305, 0, 0, 0, 0, -0.013127348731989886, 0,
      0, 0, -0.83646315449674069, 0, 0, 0, 0, 6.7819650765550188, 0, 0,
      0, 0, -1.4039352021168525, 0, 0, 0, 0, -0.2999647102418202, 0, 0, 0,
      0, 0, 0.29996471024182048, 0, 0, 0, 0, 1.403935202116853, 0, 0, 0, 0};
  // Component 3 (z), 78 elements
  const double W_ref_z[] = {
      0,
      0, 0,
      0.04151583448537164, 0.018843269537484185, 0,
      0, 0, -7.6393474522098979, 0,
      0, 0, 0.013127348731989879, 0, 0,
      0, 0, 0.0011163986474808916, 0, 0, 0,
      0, 0, -0.00064583414690861773, 0, 0, 0, 0,
      0.0011163986474808918, -0.00064583414690861859, 0, 1.0367677328022835, -0.0033118181997882339, -0.041515834485371633, -0.018843269537484185, 0,
      0, 0, -1.0367677328022835, 0, 0, 0, 0, -7.6393474522098979, 0,
      0, 0, 0.0033118181997882279, 0, 0, 0, 0, 0.013127348731989872, 0, 0,
      0, 0, 1.1307212934920607, 0, 0, 0, 0, 0.14863087674797437, 0, 0, 0,
      0, 0, -0.14863087674797437, 0, 0, 0, 0, -1.1307212934920605, 0, 0, 0, 0};
  // clang-format on

  const auto ntri = n * (n + 1) / 2;
  const double* W_refs[4] = {W_ref_s, W_ref_x, W_ref_y, W_ref_z};
  const char* comp_labels[4] = {"s", "x", "y", "z"};
  for (int comp = 0; comp < 4; ++comp) {
    auto W_comp = compute_nuclear_ltri_comp(q_engine, obs, comp);
    std::vector<double> W_ref(W_refs[comp], W_refs[comp] + ntri);
    compare_ltri(W_comp, W_ref, n, 1e-10,
                 std::string("op_q_gau_op_sap_") + comp_labels[comp]);
  }
#endif
}

TEST_CASE_METHOD(libint2::unit::DefaultFixture,
                 "op_q_gau_op point nuclear matches Operator::opVop",
                 "[engine][1-body]") {
#if defined(LIBINT2_SUPPORT_ONEBODY)
  BasisSet obs("sto-3g", atoms);
  auto point_charges = make_point_charges(atoms);
  const auto n = libint2::nbf(obs);
  const auto nshells = obs.size();
  auto shell2bf = obs.shell2bf();

  const auto lmax =
      std::min(static_cast<int>(obs.max_l()), LIBINT2_MAX_AM_elecpot);

  // Reference: Operator::opVop
  Engine ref_engine(Operator::opVop, obs.max_nprim(), lmax);
  ref_engine.set_params(point_charges);

  // op_q_gau_op with point nuclear model
  auto q_gau_data =
      libint2::make_q_gau_data(libint2::NuclearModel::PointCharge, atoms);
  Engine q_engine(Operator::op_q_gau_op,
                  std::max(obs.max_nprim(), max_nprim(q_gau_data)), lmax);
  q_engine.set_params(std::make_tuple(q_gau_data, point_charges));

  // Compare all 4 Pauli components
  for (size_t s1 = 0; s1 < nshells; ++s1) {
    auto bf1 = shell2bf[s1];
    auto n1 = obs[s1].size();
    for (size_t s2 = 0; s2 < nshells; ++s2) {
      auto bf2 = shell2bf[s2];
      auto n2 = obs[s2].size();

      ref_engine.compute(obs[s1], obs[s2]);
      const auto& ref_buf = ref_engine.results();
      q_engine.compute(obs[s1], obs[s2]);
      const auto& q_buf = q_engine.results();

      for (int comp = 0; comp < 4; ++comp) {
        if (ref_buf[comp] == nullptr && q_buf[comp] == nullptr) continue;
        REQUIRE(ref_buf[comp] != nullptr);
        REQUIRE(q_buf[comp] != nullptr);
        for (size_t i = 0; i < n1; ++i) {
          for (size_t j = 0; j < n2; ++j) {
            const auto idx = i * n2 + j;
            INFO("shell(" << s1 << "," << s2 << ") comp=" << comp << " bf("
                          << bf1 + i << "," << bf2 + j << ")");
            if (std::abs(ref_buf[comp][idx]) > 1e-12) {
              REQUIRE(q_buf[comp][idx] ==
                      Approx(ref_buf[comp][idx]).epsilon(1e-14));
            } else {
              REQUIRE(std::abs(q_buf[comp][idx]) < 1e-14);
            }
          }
        }
      }
    }
  }
#endif
}

TEST_CASE("infinite_exponent sentinel", "[engine][1-body][validation]") {
  REQUIRE(libint2::infinite_exponent ==
          std::numeric_limits<double>::infinity());
  REQUIRE(std::isinf(libint2::infinite_exponent));
  REQUIRE(libint2::infinite_exponent > 0.0);
}

TEST_CASE("make_q_gau_data factory validation",
          "[engine][1-body][validation]") {
  using libint2::Atom;

  std::vector<Atom> atoms = {Atom{1, 0.0, 0.0, 0.0}};

  SECTION("make_q_gau_data_erf rejects NaN omega") {
    REQUIRE_THROWS_AS(libint2::make_q_gau_data_erf(
                          std::numeric_limits<double>::quiet_NaN(), atoms),
                      std::invalid_argument);
  }

  SECTION("make_q_gau_data_erf rejects non-positive or non-finite omega") {
    REQUIRE_THROWS_AS(libint2::make_q_gau_data_erf(0.0, atoms),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(libint2::make_q_gau_data_erf(-1.0, atoms),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(libint2::make_q_gau_data_erf(
                          std::numeric_limits<double>::infinity(), atoms),
                      std::invalid_argument);
  }

  SECTION("make_q_gau_data_erfc rejects invalid omega") {
    REQUIRE_THROWS_AS(libint2::make_q_gau_data_erfc(
                          std::numeric_limits<double>::quiet_NaN(), atoms),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(libint2::make_q_gau_data_erfc(0.0, atoms),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(libint2::make_q_gau_data_erfc(-1.0, atoms),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(libint2::make_q_gau_data_erfc(
                          std::numeric_limits<double>::infinity(), atoms),
                      std::invalid_argument);
  }

  SECTION("make_q_gau_data_erfx rejects invalid parameters") {
    REQUIRE_THROWS_AS(
        libint2::make_q_gau_data_erfx(std::numeric_limits<double>::quiet_NaN(),
                                      0.3, 0.7, atoms),
        std::invalid_argument);
    REQUIRE_THROWS_AS(libint2::make_q_gau_data_erfx(0.0, 0.3, 0.7, atoms),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(libint2::make_q_gau_data_erfx(-1.0, 0.3, 0.7, atoms),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(
        libint2::make_q_gau_data_erfx(std::numeric_limits<double>::infinity(),
                                      0.3, 0.7, atoms),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        libint2::make_q_gau_data_erfx(
            0.5, std::numeric_limits<double>::quiet_NaN(), 0.7, atoms),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        libint2::make_q_gau_data_erfx(
            0.5, std::numeric_limits<double>::infinity(), 0.7, atoms),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        libint2::make_q_gau_data_erfx(
            0.5, 0.3, std::numeric_limits<double>::quiet_NaN(), atoms),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        libint2::make_q_gau_data_erfx(
            0.5, 0.3, std::numeric_limits<double>::infinity(), atoms),
        std::invalid_argument);
  }

  SECTION("valid factory calls succeed") {
    REQUIRE_NOTHROW(libint2::make_q_gau_data_erf(0.5, atoms));
    REQUIRE_NOTHROW(libint2::make_q_gau_data_erfc(0.5, atoms));
    REQUIRE_NOTHROW(libint2::make_q_gau_data_erfx(0.5, 0.3, 0.7, atoms));
  }
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
