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

// Tests of the zero-field London-orbital (GIAO) derivatives of the 4-center
// electron repulsion integrals (Operator::giao_coulomb) and of their geometric
// derivatives

#include <algorithm>
#include <functional>

#include "catch.hpp"
#include "fixture.h"

// see LIBINT2_ENABLE_GIAO_ERI and LIBINT2_ENABLE_GIAO_ERI3
#if defined(LIBINT2_SUPPORT_GIAO_ERI) || defined(LIBINT2_SUPPORT_GIAO_ERI3)

namespace {

using libint2::BasisSet;
using libint2::Engine;
using libint2::Operator;
using libint2::Shell;

/// all components of the integrals, [comp][element]
using Matrices = std::vector<std::vector<double>>;

/// \f$ (\vec{A} - \vec{B}) \f$ for each pair of basis functions
std::vector<std::array<double, 3>> center_differences(
    const std::vector<Shell>& shells) {
  const auto n = libint2::nbf(shells);
  std::vector<std::array<double, 3>> AB(n * n);
  size_t i = 0;
  for (const auto& s1 : shells)
    for (size_t f1 = 0; f1 != s1.size(); ++f1, ++i) {
      size_t j = 0;
      for (const auto& s2 : shells)
        for (size_t f2 = 0; f2 != s2.size(); ++f2, ++j)
          for (int xyz = 0; xyz != 3; ++xyz)
            AB[i * n + j][xyz] = s1.O[xyz] - s2.O[xyz];
    }
  return AB;
}

/// requires \p A and \p B to agree elementwise to within \p tol (absolute)
/// or 1e-10 (relative); reports the worst element only
void require_equal(const Matrices& A, const Matrices& B, double tol,
                   const std::string& label) {
  REQUIRE(A.size() == B.size());
  double worst_excess = 0;  // worst |A-B| - max(tol, 1e-10 |B|)
  size_t worst_c = 0, worst_ij = 0;
  for (size_t c = 0; c != A.size(); ++c) {
    REQUIRE(A[c].size() == B[c].size());
    for (size_t ij = 0; ij != A[c].size(); ++ij) {
      const auto excess = std::abs(A[c][ij] - B[c][ij]) -
                          std::max(tol, 1e-10 * std::abs(B[c][ij]));
      if (excess > worst_excess) {
        worst_excess = excess;
        worst_c = c;
        worst_ij = ij;
      }
    }
  }
  INFO(label << ": component " << worst_c << ", element " << worst_ij << ": "
             << A[worst_c][worst_ij] << " vs " << B[worst_c][worst_ij]);
  REQUIRE(worst_excess <= 0);
}
/// off-axis 3-center system with Cartesian s, p, d shells used to compare to
/// libcint
const std::vector<Shell>& reference_shells() {
  static const std::vector<Shell> shells = {
      Shell{{5.0, 1.2}, {{0, false, {0.4, 0.7}}}, {{0.1, -0.2, 0.05}}},
      Shell{{3.0, 0.8}, {{1, false, {0.5, 0.6}}}, {{0.1, -0.2, 0.05}}},
      Shell{{1.1}, {{2, false, {1.0}}}, {{0.1, -0.2, 0.05}}},
      Shell{{1.3, 0.3}, {{0, false, {0.5, 0.6}}}, {{1.5, 0.3, 1.1}}},
      Shell{{0.9}, {{1, false, {1.0}}}, {{1.5, 0.3, 1.1}}},
      Shell{{0.7}, {{0, false, {1.0}}}, {{-1.2, 1.3, -0.6}}},
      Shell{{0.6}, {{2, false, {1.0}}}, {{-1.2, 1.3, -0.6}}}};
  return shells;
}

/// computes all components of the 2-body integrals of \p engine over all
/// quartets of \p shells, [comp][((i * n + j) * n + k) * n + l]
Matrices compute4(Engine& engine, const std::vector<Shell>& shells) {
  const auto shell2bf = BasisSet::compute_shell2bf(shells);
  const auto n = libint2::nbf(shells);
  const auto& res = engine.results();
  Matrices M(res.size(), std::vector<double>(n * n * n * n, 0.0));
  for (size_t s1 = 0; s1 != shells.size(); ++s1)
    for (size_t s2 = 0; s2 != shells.size(); ++s2)
      for (size_t s3 = 0; s3 != shells.size(); ++s3)
        for (size_t s4 = 0; s4 != shells.size(); ++s4) {
          engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
          if (res[0] == nullptr) continue;  // screened out
          const auto n1 = shells[s1].size(), n2 = shells[s2].size(),
                     n3 = shells[s3].size(), n4 = shells[s4].size();
          for (size_t c = 0; c != res.size(); ++c)
            for (size_t i = 0, ijkl = 0; i != n1; ++i)
              for (size_t j = 0; j != n2; ++j)
                for (size_t k = 0; k != n3; ++k)
                  for (size_t l = 0; l != n4; ++l, ++ijkl)
                    M[c][((((shell2bf[s1] + i) * n + shell2bf[s2] + j) * n +
                           shell2bf[s3] + k) *
                              n +
                          shell2bf[s4] + l)] = res[c][ijkl];
        }
  return M;
}

/// the integrals of \p engine over a shell quartet, [shell set][element]
Matrices compute_quartet(Engine& engine, const std::array<Shell, 4>& s) {
  engine.compute(s[0], s[1], s[2], s[3]);
  const auto n = s[0].size() * s[1].size() * s[2].size() * s[3].size();
  Matrices result;
  for (const auto* r : engine.results()) {
    REQUIRE(r != nullptr);
    result.emplace_back(r, r + n);
  }
  return result;
}

/// requires \p computed to agree with the central finite difference
/// (\p plus - \p minus) / (2 \p h) to within \p margin (absolute) or 1e-6
/// (relative); reports the worst element only
void require_close_to_fd(const std::vector<double>& computed,
                         const std::vector<double>& plus,
                         const std::vector<double>& minus, double h,
                         double margin = 1e-8) {
  double worst_excess = 0;
  size_t worst = 0;
  for (size_t i = 0; i != computed.size(); ++i) {
    const auto ref = (plus[i] - minus[i]) / (2 * h);
    const auto excess =
        std::abs(computed[i] - ref) - std::max(margin, 1e-6 * std::abs(ref));
    if (excess > worst_excess) {
      worst_excess = excess;
      worst = i;
    }
  }
  INFO("element " << worst << ": " << computed[worst] << " vs "
                  << (plus[worst] - minus[worst]) / (2 * h));
  REQUIRE(worst_excess <= 0);
}

/// @return \p s displaced by \p d along \p xyz
Shell moved(Shell s, int xyz, double d) {
  auto O = s.O;
  O[xyz] += d;
  s.move(O);
  return s;
}

/// contracted and solid-harmonic shells on 4 centers, used to test the
/// geometric derivatives
const std::vector<Shell>& deriv_test_shells() {
  static const std::vector<Shell> shells = {
      Shell{{1.3, 0.4}, {{0, false, {0.5, 0.6}}}, {{0.1, -0.2, 0.05}}},
      Shell{{0.9, 0.35}, {{1, false, {0.7, 0.4}}}, {{1.5, 0.3, 1.1}}},
      Shell{{0.7}, {{2, true, {1.0}}}, {{-1.2, 1.3, -0.6}}},
      Shell{{0.8}, {{1, false, {1.0}}}, {{0.4, -1.1, 0.9}}}};
  return shells;
}

/// (R x e_x)_k
std::array<double, 3> cross_with_unit_vector(const std::array<double, 3>& R,
                                             size_t x) {
  std::array<double, 3> result = {{0, 0, 0}};
  result[(x + 1) % 3] = R[(x + 2) % 3];
  result[(x + 2) % 3] = -R[(x + 1) % 3];
  return result;
}

}  // namespace

#if defined(LIBINT2_SUPPORT_GIAO_ERI)

TEST_CASE("GIAO electron repulsion integrals match libcint",
          "[engine][2-body][giao]") {
  const auto& shells = reference_shells();

  // reference values computed with libcint 6.2, converted to the libint
  // conventions and normalization:
  //   X(1)(ab|cd) = -int2e_ig1(ab|cd), X(2)(ab|cd) = -int2e_ig1(cd|ab)
  // (these are also the conventions of the PySCF NMR module)
  // the quartets exercise all shell permutations done by the Engine:
  // none, bra and ket function swaps, bra-ket swap
  struct Reference {
    std::array<size_t, 4> s;     // shell indices
    std::vector<double> values;  // [comp][i][j][k][l]
  };
  const std::vector<Reference> refs = {
      {{{1, 0, 4, 3}}, {-0.000000000000000e+00, -0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        0.000000000000000e+00,  -0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        -0.000000000000000e+00, 0.000000000000000e+00,
                        0.000000000000000e+00,  -0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        0.000000000000000e+00,  0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        -0.000000000000000e+00, 0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        0.000000000000000e+00,  -0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        0.000000000000000e+00,  -0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        0.000000000000000e+00,  -0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        -0.000000000000000e+00, -0.000000000000000e+00,
                        0.000000000000000e+00,  -0.000000000000000e+00}},
      {{{0, 4, 5, 1}},
       {-2.065954639780780e-03, 1.748107884968873e-03,  -5.779538547816009e-04,
        -3.004752832249725e-03, 3.784505391540075e-03,  -1.326620812183039e-03,
        -3.712142151877517e-04, -2.070666028592435e-04, 1.859900113005574e-05,
        2.891816662586581e-03,  -2.900133009039144e-03, -1.067742973781234e-04,
        8.724632918739669e-05,  1.049526578608104e-04,  -2.770093994813214e-04,
        -2.745570390948598e-03, 4.486663821807764e-03,  -2.336566552615843e-03,
        1.377550632761716e-03,  -9.497947947017615e-04, 8.214500432221932e-04,
        3.964791238624681e-03,  -5.095984644844295e-03, 1.900736987425634e-03,
        1.802366758797287e-03,  -1.860417778000897e-03, 1.087852071167470e-03,
        5.081627480860406e-05,  -4.704384959912405e-03, -1.011515287968272e-02,
        -1.451973201262578e-04, -1.010581110967298e-03, -3.073322118072156e-03,
        -2.707626785150762e-04, -2.976859204637611e-03, -7.479382737047251e-03,
        4.602984449351487e-03,  -4.928982936260262e-04, -8.614888412861848e-03,
        1.334656529797825e-03,  -2.904258233299913e-05, -2.593053645962293e-03,
        3.105690721371375e-03,  2.222013558388264e-04,  -6.443723508243654e-03,
        1.052063925657853e-02,  8.271312319149370e-03,  3.497940373765485e-04,
        3.370371247478266e-03,  1.954140878089213e-03,  1.626742839236361e-04,
        7.708503944810248e-03,  6.466490768903284e-03,  8.863430122453294e-05}},
      {{{1, 4, 0, 3}}, {4.710258148921315e-03,  1.347407271372066e-02,
                        1.460091585942131e-03,  -1.061092429024923e-02,
                        -3.919680458138358e-03, -9.729744651068968e-03,
                        1.292913777830875e-02,  1.243009102775367e-02,
                        2.126272673334751e-03,  3.953079956810819e-03,
                        5.336876676684140e-03,  2.509640425875196e-02,
                        -3.264907967590822e-03, -4.108272944356591e-04,
                        5.301856545353345e-03,  -2.508067283285199e-02,
                        -6.167189647013270e-03, -3.107828140059753e-03,
                        -8.162763225614522e-03, -2.050680013100096e-02,
                        -1.389745747589996e-02, 1.570261713347079e-02,
                        5.421872655820505e-03,  1.044829927506656e-02,
                        -5.295672831624999e-03, -1.363669772890333e-02,
                        -1.355112069179786e-03, 4.710258148921313e-03,
                        4.771279749388771e-03,  5.604278711814460e-03,
                        3.464264492530932e-04,  -3.114481252720347e-03,
                        -3.166476577353912e-05, 7.711351711879068e-03,
                        4.039501280779865e-03,  1.321073467916741e-03,
                        1.698522181640389e-03,  1.423537126861882e-03,
                        5.274667251682687e-03,  -1.567676231865755e-04,
                        -4.108272944356587e-04, 1.157669419481015e-03,
                        -3.943777123599978e-03, -9.494035805835935e-04,
                        -8.532703648893272e-04, -7.089164285057173e-03,
                        -7.039581154833542e-03, -9.984117926077700e-03,
                        -3.872506832009930e-04, 4.348273715263157e-03,
                        -5.090514644357644e-04, -8.403813176029241e-03,
                        -4.933904764571442e-03, -1.355112069179784e-03}},
      {{{2, 3, 5, 0}}, {-2.087681400185191e-03, 1.996066337152117e-03,
                        -1.715855584329187e-03, 1.550928393364600e-03,
                        8.959297546302851e-04,  -3.524481105488151e-03,
                        -4.913814898192811e-03, -1.531244390940009e-03,
                        1.173127372182342e-03,  -1.418948718667479e-04,
                        1.950923307859170e-03,  4.758089576636779e-03,
                        5.123487056529212e-03,  -1.932257787183771e-03,
                        1.729175363780658e-03,  -2.000335537882920e-03,
                        -2.123584105154270e-03, 2.433551199395259e-03,
                        -1.392507377629879e-04, 7.965900088086791e-05,
                        1.510328078046807e-04,  2.510549622653471e-04,
                        4.493852984378240e-04,  2.323447153609037e-04,
                        -2.238757891564961e-05, -1.320680149216821e-04,
                        1.055796288942382e-04,  1.083728716016340e-04,
                        3.336683461278232e-04,  4.142483944608074e-04,
                        2.268378318744764e-04,  -4.640903438886944e-04,
                        -5.842031816111950e-05, -2.520186823730769e-04,
                        -1.287667211960560e-04, 4.912684026492867e-04}}};

  if (LIBINT2_MAX_AM_giao_eri < 2) return;  // the references include d shells
  Engine engine(Operator::giao_coulomb, 2, 2);
  for (const auto& ref : refs) {
    engine.compute(shells[ref.s[0]], shells[ref.s[1]], shells[ref.s[2]],
                   shells[ref.s[3]]);
    size_t n1234 = 1;
    for (auto s : ref.s) n1234 *= shells[s].size();
    REQUIRE(engine.results().size() == 6);
    REQUIRE(ref.values.size() == 6 * n1234);
    for (size_t c = 0; c != 6; ++c) {
      for (size_t ijkl = 0; ijkl != n1234; ++ijkl) {
        INFO("shells (" << ref.s[0] << "," << ref.s[1] << "|" << ref.s[2] << ","
                        << ref.s[3] << "), component " << c << ", element "
                        << ijkl);
        REQUIRE(
            engine.results()[c][ijkl] ==
            Approx(ref.values[c * n1234 + ijkl]).epsilon(1e-10).margin(1e-12));
      }
    }
  }
}

TEST_CASE_METHOD(libint2::unit::DefaultFixture,
                 "GIAO electron repulsion integral identities",
                 "[engine][2-body][giao]") {
  if (obs.max_l() > LIBINT2_MAX_AM_giao_eri) return;
  const std::vector<Shell>& shells = obs.shells();
  const auto n = libint2::nbf(shells);
  auto idx = [n](size_t i, size_t j, size_t k, size_t l) {
    return ((i * n + j) * n + k) * n + l;
  };
  Engine engine(Operator::giao_coulomb, obs.max_nprim(), obs.max_l());
  const auto X = compute4(engine, shells);
  REQUIRE(X.size() == 6);

  SECTION("permutational symmetry") {
    // X(1) is antisymmetric w.r.t. a<->b, symmetric w.r.t. c<->d,
    // vice versa for X(2); X(1)(ab|cd) = X(2)(cd|ab)
    auto X_ab = X, X_cd = X, X_braket = X;  // X with indices permuted
    for (size_t c = 0; c != 6; ++c) {
      const auto sign_ab = c < 3 ? -1 : 1;
      const auto sign_cd = c < 3 ? 1 : -1;
      for (size_t i = 0; i != n; ++i)
        for (size_t j = 0; j != n; ++j)
          for (size_t k = 0; k != n; ++k)
            for (size_t l = 0; l != n; ++l) {
              X_ab[c][idx(i, j, k, l)] = sign_ab * X[c][idx(j, i, k, l)];
              X_cd[c][idx(i, j, k, l)] = sign_cd * X[c][idx(i, j, l, k)];
              X_braket[c][idx(i, j, k, l)] = X[(c + 3) % 6][idx(k, l, i, j)];
            }
    }
    require_equal(X, X_ab, 1e-12, "a<->b");
    require_equal(X, X_cd, 1e-12, "c<->d");
    require_equal(X, X_braket, 1e-12, "(ab|<->|cd)");
  }

  SECTION("translation shifts GIAO ERIs by 1/2 (R x d) (ab|cd)") {
    const std::array<double, 3> d = {{0.7, -1.1, 0.4}};
    auto shifted_shells = shells;
    for (auto& s : shifted_shells)
      s.move({{s.O[0] + d[0], s.O[1] + d[1], s.O[2] + d[2]}});
    Engine eri_engine(Operator::coulomb, obs.max_nprim(), obs.max_l());
    const auto eri = compute4(eri_engine, shells)[0];
    const auto X_shifted = compute4(engine, shifted_shells);
    const auto AB = center_differences(shells);
    auto X_ref = X;
    for (size_t i = 0; i != n; ++i)
      for (size_t j = 0; j != n; ++j)
        for (size_t k = 0; k != n; ++k)
          for (size_t l = 0; l != n; ++l)
            for (int x = 0; x != 3; ++x) {
              const auto y = (x + 1) % 3, z = (x + 2) % 3;
              const auto& Rab = AB[i * n + j];
              const auto& Rcd = AB[k * n + l];
              const auto ijkl = idx(i, j, k, l);
              X_ref[x][ijkl] +=
                  0.5 * (Rab[y] * d[z] - Rab[z] * d[y]) * eri[ijkl];
              X_ref[3 + x][ijkl] +=
                  0.5 * (Rcd[y] * d[z] - Rcd[z] * d[y]) * eri[ijkl];
            }
    require_equal(X_shifted, X_ref, 1e-11, "translated giao_coulomb");
  }

  SECTION("unsupported uses throw") {
    // 2-center integrals do not depend on the magnetic field
    Engine e(Operator::giao_coulomb, 1, 1);
    REQUIRE_THROWS_AS(e.set(libint2::BraKet::xs_xs), std::invalid_argument);
    // the rejected set() leaves no kernels (rather than those of xx_xx)
    REQUIRE_THROWS_AS(e.compute(shells[0], shells[0]), Engine::missing_kernel);
#if !defined(LIBINT2_SUPPORT_GIAO_ERI3)
    Engine e3(Operator::giao_coulomb, 1, 1);
    e3.set(libint2::BraKet::xs_xx);
    REQUIRE_THROWS_AS(e3.compute(shells[0], shells[0], shells[0]),
                      Engine::missing_kernel);
#endif
#if LIBINT2_MAX_DERIV_ORDER > LIBINT2_DERIV_GIAO_ERI_ORDER
    // geometric derivatives beyond LIBINT2_ENABLE_GIAO_ERI are not generated
    Engine e1(Operator::giao_coulomb, 1, 1, LIBINT2_DERIV_GIAO_ERI_ORDER + 1);
    REQUIRE_THROWS_AS(e1.compute(shells[0], shells[0], shells[0], shells[0]),
                      Engine::missing_kernel);
#endif
  }
}

#if LIBINT2_DERIV_GIAO_ERI_ORDER > 0
TEST_CASE("geometric derivatives of GIAO electron repulsion integrals",
          "[engine][2-body][giao][deriv]") {
  // the derivatives include those of the London phase factors, check them
  // against central finite differences of the integrals;
  // contracted and solid-harmonic shells on 4 centers
  const auto& all_shells = deriv_test_shells();
  const double h = 1e-4;  // finite-difference step

  // the integrals with the center of function c / 3 displaced along
  // coordinate c % 3; N.B. only that function is displaced, even if other
  // functions of the quartet are on the same center
  auto displaced_ints = [&](Engine& engine, std::array<Shell, 4> q, size_t c,
                            double d) {
    q[c / 3] = moved(q[c / 3], c % 3, d);
    return compute_quartet(engine, q);
  };
  // all quartets of the shells with angular momenta up to lmax; these involve
  // every shell permutation done by the Engine
  auto for_each_quartet = [&](int lmax, const std::function<void(
                                            const std::array<size_t, 4>&,
                                            const std::array<Shell, 4>&)>& f) {
    std::vector<Shell> shells;
    for (const auto& s : all_shells)
      if (s.contr[0].l <= lmax) shells.push_back(s);
    const auto n = shells.size();
    for (size_t i = 0; i != n; ++i)
      for (size_t j = 0; j != n; ++j)
        for (size_t k = 0; k != n; ++k)
          for (size_t l = 0; l != n; ++l)
            f({{i, j, k, l}}, {{shells[i], shells[j], shells[k], shells[l]}});
  };

  SECTION("first derivatives") {
    const int lmax = std::min(2, static_cast<int>(LIBINT2_MAX_AM_giao_eri1));
    Engine engine0(Operator::giao_coulomb, 2, lmax, 0);
    Engine engine1(Operator::giao_coulomb, 2, lmax, 1);
    Engine eri_engine(Operator::coulomb, 2, lmax, 0);
    for_each_quartet(lmax, [&](const std::array<size_t, 4>& idx,
                               const std::array<Shell, 4>& q) {
      // shell sets are ordered as coordinate * 6 + component
      const auto d1 = compute_quartet(engine1, q);
      REQUIRE(d1.size() == 12 * 6);
      for (size_t c = 0; c != 12; ++c) {
        const auto plus = displaced_ints(engine0, q, c, h);
        const auto minus = displaced_ints(engine0, q, c, -h);
        for (size_t k = 0; k != 6; ++k) {
          INFO("shells (" << idx[0] << "," << idx[1] << "|" << idx[2] << ","
                          << idx[3] << "), coordinate " << c << ", component "
                          << k);
          require_close_to_fd(d1[c * 6 + k], plus[k], minus[k], h);
        }
      }

      // the integrals are not translationally invariant: the phase factors
      // refer to the coordinate origin, hence
      // sum_centers d X(1)_k / d x = 1/2 (R_AB x e_x)_k (ab|cd), etc.
      const auto eri = compute_quartet(eri_engine, q)[0];
      for (size_t pair = 0; pair != 2; ++pair) {
        const auto& s1 = q[2 * pair];
        const auto& s2 = q[2 * pair + 1];
        const std::array<double, 3> R = {
            {s1.O[0] - s2.O[0], s1.O[1] - s2.O[1], s1.O[2] - s2.O[2]}};
        for (size_t x = 0; x != 3; ++x) {
          const auto R_x_e = cross_with_unit_vector(R, x);
          for (size_t k = 0; k != 3; ++k) {
            for (size_t e = 0; e != eri.size(); ++e) {
              double sum = 0;
              for (size_t center = 0; center != 4; ++center)
                sum += d1[(3 * center + x) * 6 + 3 * pair + k][e];
              INFO("shells (" << idx[0] << "," << idx[1] << "|" << idx[2] << ","
                              << idx[3] << "), pair " << pair << ", x " << x
                              << ", component " << k << ", element " << e);
              REQUIRE(
                  sum ==
                  Approx(0.5 * R_x_e[k] * eri[e]).epsilon(1e-10).margin(1e-11));
            }
          }
        }
      }
    });
  }

#if LIBINT2_DERIV_GIAO_ERI_ORDER > 1
  SECTION("second derivatives") {
    const int lmax = std::min(2, static_cast<int>(LIBINT2_MAX_AM_giao_eri2));
    Engine engine1(Operator::giao_coulomb, 2, lmax, 1);
    Engine engine2(Operator::giao_coulomb, 2, lmax, 2);
    for_each_quartet(lmax, [&](const std::array<size_t, 4>& idx,
                               const std::array<Shell, 4>& q) {
      // upper triangle of the 12x12 coordinates, ordered as
      // pair * 6 + component
      const auto d2 = compute_quartet(engine2, q);
      REQUIRE(d2.size() == 78 * 6);
      for (size_t c2 = 0; c2 != 12; ++c2) {
        const auto plus = displaced_ints(engine1, q, c2, h);
        const auto minus = displaced_ints(engine1, q, c2, -h);
        for (size_t c1 = 0; c1 != 12; ++c1) {
          const auto i = std::min(c1, c2), j = std::max(c1, c2);
          const auto c12 = i * (2 * 12 - i - 1) / 2 + j;
          for (size_t k = 0; k != 6; ++k) {
            INFO("shells (" << idx[0] << "," << idx[1] << "|" << idx[2] << ","
                            << idx[3] << "), coordinates (" << c1 << "," << c2
                            << "), component " << k);
            // N.B. looser: the differentiated functions vary faster
            require_close_to_fd(d2[c12 * 6 + k], plus[c1 * 6 + k],
                                minus[c1 * 6 + k], h, 1e-7);
          }
        }
      }
    });
  }
#endif  // LIBINT2_DERIV_GIAO_ERI_ORDER > 1
}
#endif  // LIBINT2_DERIV_GIAO_ERI_ORDER > 0

#endif  // defined(LIBINT2_SUPPORT_GIAO_ERI)

#if defined(LIBINT2_SUPPORT_GIAO_ERI3)

namespace {

/// the integrals of \p engine over a shell triplet, [shell set][element]
Matrices compute_triplet(Engine& engine, const std::array<Shell, 3>& s,
                         libint2::BraKet braket = libint2::BraKet::xs_xx) {
  if (braket == libint2::BraKet::xs_xx)
    engine.compute(s[0], s[1], s[2]);
  else  // xx_xs: the unpaired function is last
    engine.compute(s[1], s[2], s[0]);
  const auto n = s[0].size() * s[1].size() * s[2].size();
  Matrices result;
  for (const auto* r : engine.results()) {
    REQUIRE(r != nullptr);
    result.emplace_back(r, r + n);
  }
  return result;
}

/// calls \p f for all triplets of the test shells with angular momenta up to
/// \p lmax ; these involve every shell permutation done by the Engine
void for_each_triplet(
    int lmax, const std::function<void(const std::array<size_t, 3>&,
                                       const std::array<Shell, 3>&)>& f) {
  std::vector<Shell> shells;
  for (const auto& s : deriv_test_shells())
    if (s.contr[0].l <= lmax) shells.push_back(s);
  const auto n = shells.size();
  for (size_t i = 0; i != n; ++i)
    for (size_t j = 0; j != n; ++j)
      for (size_t k = 0; k != n; ++k)
        f({{i, j, k}}, {{shells[i], shells[j], shells[k]}});
}

/// an Engine for the 3-center integrals of \p oper
Engine make_engine3(Operator oper, int lmax, int deriv_order,
                    libint2::BraKet braket = libint2::BraKet::xs_xx) {
  return Engine(
      oper, 2, lmax, deriv_order, std::numeric_limits<double>::epsilon(),
      libint2::operator_traits<Operator::coulomb>::default_params(), braket);
}

}  // namespace

TEST_CASE("3-center GIAO electron repulsion integrals",
          "[engine][2-body][giao]") {
  // the fitting function X has no London phase factor, hence
  // X_k = 1/2 (R_CD x (X|c r d))_k
  const int lmax =
      std::min(2, std::min(static_cast<int>(LIBINT2_MAX_AM_3giao_eri),
                           static_cast<int>(LIBINT2_MAX_AM_default)));
  auto engine = make_engine3(Operator::giao_coulomb, lmax, 0);
  auto engine_xx_xs =
      make_engine3(Operator::giao_coulomb, lmax, 0, libint2::BraKet::xx_xs);
  auto eri_engine = make_engine3(Operator::coulomb, lmax, 0);

  for_each_triplet(lmax, [&](const std::array<size_t, 3>& idx,
                             const std::array<Shell, 3>& t) {
    INFO("shells (" << idx[0] << "|" << idx[1] << "," << idx[2] << ")");
    const auto X = compute_triplet(engine, t);
    REQUIRE(X.size() == 3);
    const auto n0 = t[0].size(), n1 = t[1].size(), n2 = t[2].size();

    // antisymmetric w.r.t. c<->d
    const auto X_dc = compute_triplet(engine, {{t[0], t[2], t[1]}});
    // BraKet::xx_xs: (cd|X)
    const auto X_xx_xs =
        compute_triplet(engine_xx_xs, t, libint2::BraKet::xx_xs);
    REQUIRE(X_xx_xs.size() == 3);
    for (size_t k = 0; k != 3; ++k)
      for (size_t x = 0; x != n0; ++x)
        for (size_t c = 0; c != n1; ++c)
          for (size_t d = 0; d != n2; ++d) {
            const auto value = X[k][(x * n1 + c) * n2 + d];
            REQUIRE(value == Approx(-X_dc[k][(x * n2 + d) * n1 + c])
                                 .epsilon(1e-10)
                                 .margin(1e-12));
            REQUIRE(value == Approx(X_xx_xs[k][(c * n2 + d) * n0 + x])
                                 .epsilon(1e-10)
                                 .margin(1e-12));
          }

    // translation by d shifts X by 1/2 (R_CD x d) (X|cd)
    const std::array<double, 3> d = {{0.7, -1.1, 0.4}};
    auto t_shifted = t;
    for (auto& s : t_shifted)
      for (int xyz = 0; xyz != 3; ++xyz) s = moved(s, xyz, d[xyz]);
    const auto X_shifted = compute_triplet(engine, t_shifted);
    const auto eri = compute_triplet(eri_engine, t)[0];
    const std::array<double, 3> R = {
        {t[1].O[0] - t[2].O[0], t[1].O[1] - t[2].O[1], t[1].O[2] - t[2].O[2]}};
    for (size_t k = 0; k != 3; ++k) {
      const auto l = (k + 1) % 3, m = (k + 2) % 3;
      const auto R_x_d = R[l] * d[m] - R[m] * d[l];
      for (size_t e = 0; e != eri.size(); ++e)
        REQUIRE(X_shifted[k][e] == Approx(X[k][e] + 0.5 * R_x_d * eri[e])
                                       .epsilon(1e-10)
                                       .margin(1e-11));
    }

#if defined(LIBINT2_SUPPORT_GIAO_ERI)
    // = the components of the ket pair of the 4-center integrals with a unit
    // function in the bra
    if (lmax <= static_cast<int>(LIBINT2_MAX_AM_giao_eri)) {
      Engine engine4(Operator::giao_coulomb, 2, lmax);
      const auto X4 =
          compute_quartet(engine4, {{t[0], Shell::unit(), t[1], t[2]}});
      REQUIRE(X4.size() == 6);
      for (size_t k = 0; k != 3; ++k)
        for (size_t e = 0; e != X[k].size(); ++e)
          REQUIRE(X[k][e] == Approx(X4[3 + k][e]).epsilon(1e-10).margin(1e-12));
    }
#endif
  });

  // 2-center integrals do not depend on the magnetic field
  REQUIRE_THROWS_AS(engine.set(libint2::BraKet::xs_xs), std::invalid_argument);
  // the rejected set() leaves no kernels (rather than those of xs_xx)
  const auto& s0 = deriv_test_shells()[0];
  REQUIRE_THROWS_AS(engine.compute(s0, s0), Engine::missing_kernel);
}

#if LIBINT2_DERIV_GIAO_ERI3_ORDER > 0
TEST_CASE("geometric derivatives of 3-center GIAO electron repulsion integrals",
          "[engine][2-body][giao][deriv]") {
  const double h = 1e-4;  // finite-difference step
  // the integrals with the center of function c / 3 displaced along
  // coordinate c % 3
  auto displaced_ints = [&](Engine& engine, std::array<Shell, 3> t, size_t c,
                            double d) {
    t[c / 3] = moved(t[c / 3], c % 3, d);
    return compute_triplet(engine, t);
  };

  SECTION("first derivatives") {
    const int lmax =
        std::min(2, std::min(static_cast<int>(LIBINT2_MAX_AM_3giao_eri1),
                             static_cast<int>(LIBINT2_MAX_AM_default1)));
    auto engine0 = make_engine3(Operator::giao_coulomb, lmax, 0);
    auto engine1 = make_engine3(Operator::giao_coulomb, lmax, 1);
    auto eri_engine = make_engine3(Operator::coulomb, lmax, 0);
    for_each_triplet(lmax, [&](const std::array<size_t, 3>& idx,
                               const std::array<Shell, 3>& t) {
      // shell sets are ordered as coordinate * 3 + component
      const auto d1 = compute_triplet(engine1, t);
      REQUIRE(d1.size() == 9 * 3);
      for (size_t c = 0; c != 9; ++c) {
        const auto plus = displaced_ints(engine0, t, c, h);
        const auto minus = displaced_ints(engine0, t, c, -h);
        for (size_t k = 0; k != 3; ++k) {
          INFO("shells (" << idx[0] << "|" << idx[1] << "," << idx[2]
                          << "), coordinate " << c << ", component " << k);
          require_close_to_fd(d1[c * 3 + k], plus[k], minus[k], h);
        }
      }

      // sum_centers d X_k / d x = 1/2 (R_CD x e_x)_k (X|cd)
      const auto eri = compute_triplet(eri_engine, t)[0];
      const std::array<double, 3> R = {{t[1].O[0] - t[2].O[0],
                                        t[1].O[1] - t[2].O[1],
                                        t[1].O[2] - t[2].O[2]}};
      for (size_t x = 0; x != 3; ++x) {
        const auto R_x_e = cross_with_unit_vector(R, x);
        for (size_t k = 0; k != 3; ++k)
          for (size_t e = 0; e != eri.size(); ++e) {
            double sum = 0;
            for (size_t center = 0; center != 3; ++center)
              sum += d1[(3 * center + x) * 3 + k][e];
            INFO("shells (" << idx[0] << "|" << idx[1] << "," << idx[2]
                            << "), x " << x << ", component " << k
                            << ", element " << e);
            REQUIRE(
                sum ==
                Approx(0.5 * R_x_e[k] * eri[e]).epsilon(1e-10).margin(1e-11));
          }
      }
    });
  }

#if LIBINT2_DERIV_GIAO_ERI3_ORDER > 1
  SECTION("second derivatives") {
    const int lmax =
        std::min(2, std::min(static_cast<int>(LIBINT2_MAX_AM_3giao_eri2),
                             static_cast<int>(LIBINT2_MAX_AM_default2)));
    auto engine1 = make_engine3(Operator::giao_coulomb, lmax, 1);
    auto engine2 = make_engine3(Operator::giao_coulomb, lmax, 2);
    for_each_triplet(lmax, [&](const std::array<size_t, 3>& idx,
                               const std::array<Shell, 3>& t) {
      // upper triangle of the 9x9 coordinates, ordered as pair * 3 + component
      const auto d2 = compute_triplet(engine2, t);
      REQUIRE(d2.size() == 45 * 3);
      for (size_t c2 = 0; c2 != 9; ++c2) {
        const auto plus = displaced_ints(engine1, t, c2, h);
        const auto minus = displaced_ints(engine1, t, c2, -h);
        for (size_t c1 = 0; c1 != 9; ++c1) {
          const auto i = std::min(c1, c2), j = std::max(c1, c2);
          const auto c12 = i * (2 * 9 - i - 1) / 2 + j;
          for (size_t k = 0; k != 3; ++k) {
            INFO("shells (" << idx[0] << "|" << idx[1] << "," << idx[2]
                            << "), coordinates (" << c1 << "," << c2
                            << "), component " << k);
            // N.B. looser: the differentiated functions vary faster
            require_close_to_fd(d2[c12 * 3 + k], plus[c1 * 3 + k],
                                minus[c1 * 3 + k], h, 1e-7);
          }
        }
      }
    });
  }
#endif  // LIBINT2_DERIV_GIAO_ERI3_ORDER > 1
}
#endif  // LIBINT2_DERIV_GIAO_ERI3_ORDER > 0

#endif  // defined(LIBINT2_SUPPORT_GIAO_ERI3)

#endif  // LIBINT2_SUPPORT_GIAO_ERI || LIBINT2_SUPPORT_GIAO_ERI3
