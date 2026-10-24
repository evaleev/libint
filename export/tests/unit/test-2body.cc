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

#include <libint2/deriv_iter.h>

#include <algorithm>
#include <limits>
#include <sstream>

#include "catch.hpp"
#include "fixture.h"

#if defined(NO_LIBINT_COMPILER_CODE)
#include "../eri/eri.h"
#else
#include <test_eri/eri.h>
#endif

typedef unsigned int uint;

TEST_CASE("Slater/Yukawa integrals", "[engine][2-body]") {
#if defined(LIBINT2_SUPPORT_ERI)
  std::vector<Shell> obs{
      // pseudorandom d
      Shell{{1.0, 3.0}, {{2, true, {1.0, 0.3}}}, {{0.0, 0.0, 0.0}}},
      // pseudorandom d
      Shell{{2.0, 5.0}, {{2, true, {1.0, 0.2}}}, {{1.0, 1.0, 1.0}}},
      //    tight functions are problematic with Gaussian fits. The errors in
      //    some integrals involving the functions below can be huge
      //      // O 1s in cc-pVDZ: tightest primitive
      //      Shell{{11720.0000000},
      //            {{0,
      //              false,
      //              {1.0}}},
      //            {{2.0, 2.0, -1.0}}},
      //      // O 1s in STO-3G: tightest primitive
      //      Shell{{130.709320000},
      //            {{0, false, {1.0}}},
      //            {{-1.0, -1.0, 0.0}}},
      //      // O 1s in cc-pVDZ
      //      Shell{{11720.0000000, 1759.0000000, 400.8000000,
      //      113.7000000, 37.0300000,
      //             13.2700000, 5.0250000, 1.0130000},
      //            {{0,
      //              false,
      //              {0.0007100, 0.0054700, 0.0278370, 0.1048000, 0.2830620,
      //              0.4487190,
      //               0.2709520, 0.0154580}}},
      //            {{2.0, 2.0, -1.0}}},
      //      // O 1s in STO-3G
      //      Shell{{130.709320000, 23.808861000, 6.443608300},
      //            {{0, false, {0.15432897, 0.53532814, 0.44463454}}},
      //            {{-1.0, -1.0, 0.0}}}
  };
  const auto max_nprim = libint2::max_nprim(obs);
  const auto max_l = libint2::max_l(obs);

  // 6-term GTG fit of exp(-r12) from DOI 10.1063/1.1999632
  //  std::vector<std::pair<double, double>>
  //  cgtg_params{{0.2209,0.3144},{1.004,0.3037},{3.622,0.1681},{12.16,0.09811},{45.87,0.06024},{254.4,0.03726}};
  // "precise" fit of exp(-r12) on [0,10)
  std::vector<std::pair<double, double>> cgtg_params = {
      {0.10535330565471572, 0.08616353459042002},
      {0.22084823136587992, 0.08653979627551414},
      {0.3543431104992702, 0.08803697599356214},
      {0.48305514665749105, 0.09192519612306953},
      {0.6550035700167584, 0.10079776426873248},
      {1.1960917050801643, 0.11666110644901695},
      {2.269278814810891, 0.14081371547404428},
      {5.953990617813977, 0.13410255216448014},
      {18.31911063199608, 0.0772095196191394},
      {66.98443868169818, 0.049343985939540556},
      {367.24137290439205, 0.03090625839896873},
      {5.655142311118115, -0.017659052507938647}};

  if (LIBINT2_MAX_AM_eri >= max_l) {
    for (int k = -1; k <= 0; ++k) {
      auto engine =
          Engine(k == 0 ? Operator::stg : Operator::yukawa, max_nprim, max_l);
      const auto scale = 2.3;
      engine.set_params(1.0).prescale_by(scale);
      REQUIRE(engine.prescaled_by() == scale);
      const auto &results = engine.results();
      auto engine_ref = Engine(
          k == 0 ? Operator::cgtg : Operator::cgtg_x_coulomb, max_nprim, max_l);
      engine_ref.set_params(cgtg_params);
      const auto &results_ref = engine_ref.results();

      const auto nshell = obs.size();
      for (int s1 = 0; s1 != nshell; ++s1) {
        for (int s2 = 0; s2 != nshell; ++s2) {
          for (int s3 = 0; s3 != nshell; ++s3) {
            for (int s4 = 0; s4 != nshell; ++s4) {
              engine.compute(obs[s1], obs[s2], obs[s3], obs[s4]);
              engine_ref.compute(obs[s1], obs[s2], obs[s3], obs[s4]);
              if (results[0] != nullptr) {
                REQUIRE(results_ref[0] != nullptr);
                const auto setsize = obs[s1].size() * obs[s2].size() *
                                     obs[s3].size() * obs[s4].size();
                for (int i = 0; i != setsize; ++i) {
                  REQUIRE(results[0][i] / scale ==
                          Approx(results_ref[0][i]).margin(1e-3));
                }
              }
            }
          }
        }
      }
    }
  }
#endif  // LIBINT2_SUPPORT_ERI
}

// see https://github.com/evaleev/libint/issues/216
TEST_CASE_METHOD(libint2::unit::DefaultFixture, "bra-ket permutation",
                 "[engine][2-body]") {
#if defined(LIBINT2_SUPPORT_ERI)
  if (LIBINT2_MAX_AM_eri < 2) return;
  using libint2::BasisSet;
  using libint2::Engine;
  using libint2::Operator;
  using libint2::Shell;
  std::vector<Shell> shells;
  shells.push_back({
      {0.8378385011e+02, 0.1946956493e+02, 0.6332106784e+01},  // exponents
      {// P shell, spherical=false, contraction coefficients
       {1, false, {0.1559162750, 0.6076837186, 0.3919573931}}},
      {{0, 0, 0}}  // origin coordinates
  });
  shells.push_back({
      {0.2964817927e+01, 0.9043639676e+00, 0.3489317337e+00},  // exponents
      {// D shell, spherical=false, contraction coefficients
       {2, false, {0.2197679508, 0.6555473627, 0.2865732590}}},
      {{0, 0, 0}}  // origin coordinates
  });

  auto obs = BasisSet(std::move(shells));

  auto engine = Engine(libint2::Operator::coulomb, libint2::max_nprim(obs),
                       libint2::max_l(obs), 0);
  auto engine_kb = Engine(libint2::Operator::coulomb, libint2::max_nprim(obs),
                          libint2::max_l(obs), 0);
  // Force uniform normalised Cartesian functions
  engine.set(libint2::CartesianShellNormalization::uniform);
  engine_kb.set(libint2::CartesianShellNormalization::uniform);
  const auto &buf = engine.results();
  const auto &buf_kb = engine_kb.results();

  for (auto s1 = 0; s1 != obs.size(); ++s1) {
    auto n1 = obs[s1].size();  // number of basis functions in shell 1
    for (auto s2 = 0; s2 != obs.size(); ++s2) {
      auto n2 = obs[s2].size();  // number of basis functions in shell 2
      for (auto s3 = 0; s3 != obs.size(); ++s3) {
        auto n3 = obs[s3].size();  // number of basis functions in shell 3
        for (auto s4 = 0; s4 != obs.size(); ++s4) {
          auto n4 = obs[s4].size();  // number of basis functions in shell 4
          engine.compute2<Operator::coulomb, libint2::BraKet::xx_xx, 0>(
              obs[s1], obs[s2], obs[s3], obs[s4]);
          const auto *buf_1234 = buf[0];
          engine_kb.compute2<Operator::coulomb, libint2::BraKet::xx_xx, 0>(
              obs[s3], obs[s4], obs[s1], obs[s2]);
          const auto *buf_3412 = buf_kb[0];
          for (auto f1 = 0ul, f1234 = 0ul; f1 != n1; ++f1) {
            for (auto f2 = 0ul; f2 != n2; ++f2) {
              for (auto f3 = 0ul; f3 != n3; ++f3) {
                for (auto f4 = 0ul; f4 != n4; ++f4, ++f1234) {
                  const auto integral = buf_1234[f1234];
                  const auto f3412 = ((f3 * n4 + f4) * n1 + f1) * n2 + f2;
                  const auto integral_kb = buf_3412[f3412];
                  REQUIRE(std::abs(integral - integral_kb) < 1e-14);
                }
              }
            }
          }
        }
      }
    }
  }
#endif  // LIBINT2_SUPPORT_ERI
}

TEST_CASE("eri geometric derivatives", "[engine][2-body]") {
#if defined(LIBINT2_SUPPORT_ERI)
  std::vector<Shell> obs{
      // pseudorandom p
      Shell{{1.0, 0.3}, {{1, false, {0.9, 0.3}}}, {{0.0, 0.0, 0.0}}},
      // pseudorandom p
      Shell{{2.0, 0.4}, {{1, false, {0.8, -0.2}}}, {{1.0, 1.0, 1.0}}},
  };
  const auto max_nprim = libint2::max_nprim(obs);
  const auto max_l = libint2::max_l(obs);

  {
    const auto deriv_order = LIBINT_INCLUDE_ERI;
    Engine engine;
    try {
      engine = Engine(Operator::coulomb, max_nprim, max_l, deriv_order);
    } catch (Engine::lmax_exceeded &) {  // skip the test if lmax exceeded
      return;
    }
    const auto &buf = engine.results();

    libint2::CartesianDerivIterator<4> diter(deriv_order);
    const unsigned int nderiv = diter.range_size();

    const auto nshell = obs.size();
    for (int s0 = 0; s0 != nshell; ++s0) {
      for (int s1 = 0; s1 != nshell; ++s1) {
        for (int s2 = 0; s2 != nshell; ++s2) {
          for (int s3 = 0; s3 != nshell; ++s3) {
            engine.compute(obs[s0], obs[s1], obs[s2], obs[s3]);

            {
              // compare Libint integrals against the reference method
              // since the reference implementation computes integrals one at a
              // time (not one shell-set at a time) the outer loop is over the
              // basis functions

              LIBINT2_REF_REALTYPE Aref[3];
              for (int i = 0; i < 3; ++i) Aref[i] = obs[s0].O[i];
              LIBINT2_REF_REALTYPE Bref[3];
              for (int i = 0; i < 3; ++i) Bref[i] = obs[s1].O[i];
              LIBINT2_REF_REALTYPE Cref[3];
              for (int i = 0; i < 3; ++i) Cref[i] = obs[s2].O[i];
              LIBINT2_REF_REALTYPE Dref[3];
              for (int i = 0; i < 3; ++i) Dref[i] = obs[s3].O[i];

              int ijkl = 0;

              int l0, m0, n0;
              FOR_CART(l0, m0, n0, obs[s0].contr[0].l)

              int l1, m1, n1;
              FOR_CART(l1, m1, n1, obs[s1].contr[0].l)

              int l2, m2, n2;
              FOR_CART(l2, m2, n2, obs[s2].contr[0].l)

              int l3, m3, n3;
              FOR_CART(l3, m3, n3, obs[s3].contr[0].l)

              {
                //
                // compute reference integrals
                //
                std::vector<LIBINT2_REF_REALTYPE> ref_eri(nderiv, 0.0);

                uint p0123 = 0;
                for (uint p0 = 0; p0 < obs[s0].nprim(); p0++) {
                  for (uint p1 = 0; p1 < obs[s1].nprim(); p1++) {
                    for (uint p2 = 0; p2 < obs[s2].nprim(); p2++) {
                      for (uint p3 = 0; p3 < obs[s3].nprim(); p3++, p0123++) {
                        const LIBINT2_REF_REALTYPE alpha0 = obs[s0].alpha[p0];
                        const LIBINT2_REF_REALTYPE alpha1 = obs[s1].alpha[p1];
                        const LIBINT2_REF_REALTYPE alpha2 = obs[s2].alpha[p2];
                        const LIBINT2_REF_REALTYPE alpha3 = obs[s3].alpha[p3];

                        const LIBINT2_REF_REALTYPE c0 =
                            obs[s0].contr[0].coeff[p0];
                        const LIBINT2_REF_REALTYPE c1 =
                            obs[s1].contr[0].coeff[p1];
                        const LIBINT2_REF_REALTYPE c2 =
                            obs[s2].contr[0].coeff[p2];
                        const LIBINT2_REF_REALTYPE c3 =
                            obs[s3].contr[0].coeff[p3];
                        const LIBINT2_REF_REALTYPE c0123 = c0 * c1 * c2 * c3;

                        libint2::CartesianDerivIterator<4> diter(deriv_order);
                        bool last_deriv = false;
                        unsigned int di = 0;
                        do {
                          ref_eri[di++] +=
                              c0123 * eri(&(*diter)[0], l0, m0, n0, alpha0,
                                          Aref, l1, m1, n1, alpha1, Bref, l2,
                                          m2, n2, alpha2, Cref, l3, m3, n3,
                                          alpha3, Dref, 0);
                          last_deriv = diter.last();
                          if (!last_deriv) diter.next();
                        } while (!last_deriv);
                      }
                    }
                  }
                }

                //
                // extract Libint integrals
                //
                std::vector<LIBINT2_REALTYPE> new_eri;
                for (auto d = 0; d != nderiv; ++d)
                  new_eri.push_back(buf[d][ijkl]);

                //
                // compare reference and libint integrals
                //
                for (unsigned int di = 0; di < nderiv; ++di) {
                  const double ABSOLUTE_DEVIATION_THRESHOLD =
                      5.0E-14;  // indicate failure if any integral differs in
                                // absolute sense by more than this loss of
                                // precision in HRR likely limits precision for
                                // high-L (e.g. (dp|dd), (dd|dd), etc.)
                  const double RELATIVE_DEVIATION_THRESHOLD =
                      1.0E-9;  // indicate failure if any integral differs in
                               // relative sense by more than this

                  const LIBINT2_REF_REALTYPE abs_error =
                      abs(ref_eri[di] - LIBINT2_REF_REALTYPE(new_eri[di]));
                  const LIBINT2_REF_REALTYPE relabs_error =
                      abs(abs_error / ref_eri[di]);
                  bool not_ok =
                      relabs_error > RELATIVE_DEVIATION_THRESHOLD &&
                      abs_error > ABSOLUTE_DEVIATION_THRESHOLD *
                                      std::pow(3., deriv_order > 2
                                                       ? deriv_order - 2
                                                       : 0);
                  if (not_ok) {
                    std::cout << "Elem " << ijkl << " di= " << di
                              << " : ref = " << ref_eri[di]
                              << " libint = " << new_eri[di]
                              << " relabs_error = " << relabs_error
                              << " abs_error = " << abs_error << std::endl;
                  }
                  REQUIRE(!not_ok);
                }
              }
              ++ijkl;
              END_FOR_CART
              END_FOR_CART
              END_FOR_CART
              END_FOR_CART

            }  // checking computed values vs. the reference
          }
        }
      }
    }
  }
#endif  // LIBINT2_SUPPORT_ERI
}

// Validates geometric derivatives of Coulomb integrals against 5-point central
// finite differences of the derivatives of the next lower order, for every
// combination of s, p, d shells (and up to f on the fitting center), with
// contracted and uncontracted shells, in all brakets with derivative kernels.
// Complements "eri geometric derivatives" above, which compares (pp|pp) against
// an independent reference: this covers every kernel in the build, including
// the code-sharing paths that only appear with spectator shells of unequal
// size, at the cost of relying on the correctness of the next lower order.
namespace {

// index of the second derivative with respect to coordinates i and j: second
// derivatives are the upper triangle of the Hessian, in row-major order
size_t hessian_index(size_t ncoords, size_t i, size_t j) {
  const auto row = std::min(i, j);
  const auto col = std::max(i, j);
  return row * (2 * ncoords - row - 1) / 2 + col;
}

// all centers distinct and in general position
Shell make_fd_shell(int l, bool contracted, int center) {
  const std::array<double, 3> O[] = {{{0.0, 0.0, 0.0}},
                                     {{0.7, 0.2, -0.3}},
                                     {{-0.4, 0.9, 0.5}},
                                     {{0.3, -0.6, 1.1}}};
  const double scale = 1.0 + 0.2 * center;
  if (contracted)
    return Shell{{2.3 * scale, 0.8 * scale, 0.31 * scale},
                 {{l, false, {0.5, 0.4, 0.3}}},
                 O[center]};
  return Shell{{0.65 * scale}, {{l, false, {1.0}}}, O[center]};
}

std::vector<std::vector<double>> compute_fd_shellsets(
    Engine &engine, BraKet braket, const std::vector<Shell> &shells) {
  const Engine::target_ptr_vec *results = nullptr;
  switch (braket) {
    case BraKet::xx_xx:
      results = &engine.compute(shells[0], shells[1], shells[2], shells[3]);
      break;
    case BraKet::xs_xx:
      results = &engine.compute(shells[0], shells[1], shells[2]);
      break;
    case BraKet::xs_xs:
      results = &engine.compute(shells[0], shells[1]);
      break;
    default:
      FAIL("unsupported braket");
  }
  size_t n = 1;
  for (const auto &sh : shells) n *= sh.size();
  std::vector<std::vector<double>> shellsets;
  for (const auto *target : *results) {
    if (target == nullptr)  // screened out
      shellsets.emplace_back(n, 0.0);
    else
      shellsets.emplace_back(target, target + n);
  }
  return shellsets;
}

// checks derivatives of order deriv_order of a shell set with the given
// angular momenta against finite differences of order deriv_order - 1;
// returns false if the library lacks the kernels for this shell set
bool check_fd_derivs(BraKet braket, int deriv_order, const std::vector<int> &l,
                     const std::vector<bool> &contracted) {
  const auto ncenters = l.size();
  std::vector<Shell> shells;
  for (size_t c = 0; c != ncenters; ++c)
    shells.push_back(make_fd_shell(l[c], contracted[c], c));
  const auto nprim = libint2::max_nprim(shells);
  const auto lmax = libint2::max_l(shells);
  const auto eps = std::numeric_limits<double>::epsilon();
  Engine engine, engine_lower;
  try {
    engine = Engine(Operator::coulomb, nprim, lmax, deriv_order, eps,
                    libint2::default_params(Operator::coulomb), braket);
    engine_lower = Engine(Operator::coulomb, nprim, lmax, deriv_order - 1, eps,
                          libint2::default_params(Operator::coulomb), braket);
  } catch (Engine::lmax_exceeded &) {
    return false;
  }
  engine.set_precision(0.0);
  engine_lower.set_precision(0.0);

  const auto analytic = compute_fd_shellsets(engine, braket, shells);
  const auto ncoords = 3 * ncenters;
  REQUIRE(analytic.size() ==
          libint2::num_geometrical_derivatives(ncenters, deriv_order));

  const double h = 1e-3;
  for (size_t coord1 = 0; coord1 != ncoords; ++coord1) {
    // derivatives of order deriv_order - 1 at a displaced geometry;
    // N.B. displace Shell::O in place, Shell's constructor would renormalize
    // the contraction coefficients
    auto lower = [&](double delta) {
      auto displaced = shells;
      displaced[coord1 / 3].O[coord1 % 3] += delta;
      return compute_fd_shellsets(engine_lower, braket, displaced);
    };
    const auto m2 = lower(-2 * h);
    const auto m1 = lower(-h);
    const auto p1 = lower(h);
    const auto p2 = lower(2 * h);

    for (size_t s = 0; s != p1.size(); ++s) {
      const auto deriv_index =
          deriv_order == 1 ? coord1 : hessian_index(ncoords, s, coord1);
      const auto &a = analytic.at(deriv_index);
      double max_error = 0, max_numeric = 0;
      for (size_t i = 0; i != a.size(); ++i) {
        const auto numeric =
            (m2[s][i] - 8 * m1[s][i] + 8 * p1[s][i] - p2[s][i]) / (12 * h);
        max_error = std::max(max_error, std::abs(a[i] - numeric));
        max_numeric = std::max(max_numeric, std::abs(numeric));
      }
      std::ostringstream oss;
      for (size_t c = 0; c != ncenters; ++c)
        oss << (contracted[c] ? "LPDF"[l[c]] : "lpdf"[l[c]]);
      INFO("braket=" << static_cast<int>(braket)
                     << " deriv_order=" << deriv_order
                     << " shells=" << oss.str() << " coord=" << coord1
                     << " lower_shellset=" << s << " max_error=" << max_error);
      // 5-point stencil with h=1e-3: truncation and roundoff errors are both
      // ~1e-12 relative to the scale of the derivatives
      REQUIRE(max_error <= 1e-9 * std::max(1.0, max_numeric));
    }
  }
  return true;
}

// contraction pattern for a shell set, cycling through the patterns so that
// each angular momentum combination is tested once
std::vector<bool> fd_contraction_pattern(const std::vector<int> &l) {
  int key = 0;
  for (size_t c = 0; c != l.size(); ++c) key += (c + 1) * l[c];
  std::vector<bool> contracted(l.size(), false);
  switch (key % 3) {
    case 0:
      break;
    case 1:
      for (size_t c = 0; c < l.size(); c += 2) contracted[c] = true;
      break;
    case 2:
      contracted.assign(l.size(), true);
      break;
  }
  return contracted;
}

}  // namespace

TEST_CASE("eri geometric derivatives vs finite differences",
          "[engine][2-body][deriv]") {
  size_t ntested = 0;
#if defined(LIBINT2_SUPPORT_ERI) && LIBINT2_DERIV_ERI_ORDER > 0
  SECTION("4-center") {
    const int lmax = std::min(2, LIBINT2_MAX_AM_eri);
    for (int d = 1; d <= LIBINT2_DERIV_ERI_ORDER; ++d)
      for (int la = 0; la <= lmax; ++la)
        for (int lb = 0; lb <= lmax; ++lb)
          for (int lc = 0; lc <= lmax; ++lc)
            for (int ld = 0; ld <= lmax; ++ld) {
              const std::vector<int> l{la, lb, lc, ld};
              if (check_fd_derivs(BraKet::xx_xx, d, l,
                                  fd_contraction_pattern(l)))
                ++ntested;
            }
    REQUIRE(ntested > 0);
  }
#endif
#if defined(LIBINT2_SUPPORT_ERI3) && LIBINT2_DERIV_ERI3_ORDER > 0
  SECTION("3-center") {
    const int lmax_fit = std::min(3, LIBINT2_MAX_AM_3eri);
    const int lmax = std::min(2, LIBINT2_MAX_AM_default);
    for (int d = 1; d <= LIBINT2_DERIV_ERI3_ORDER; ++d)
      for (int lp = 0; lp <= lmax_fit; ++lp)
        for (int lc = 0; lc <= lmax; ++lc)
          for (int ld = 0; ld <= lmax; ++ld) {
            const std::vector<int> l{lp, lc, ld};
            if (check_fd_derivs(BraKet::xs_xx, d, l, fd_contraction_pattern(l)))
              ++ntested;
          }
    REQUIRE(ntested > 0);
  }
#endif
#if defined(LIBINT2_SUPPORT_ERI2) && LIBINT2_DERIV_ERI2_ORDER > 0
  SECTION("2-center") {
    const int lmax = std::min(3, LIBINT2_MAX_AM_2eri);
    for (int d = 1; d <= LIBINT2_DERIV_ERI2_ORDER; ++d)
      for (int lp = 0; lp <= lmax; ++lp)
        for (int lq = 0; lq <= lmax; ++lq) {
          const std::vector<int> l{lp, lq};
          if (check_fd_derivs(BraKet::xs_xs, d, l, fd_contraction_pattern(l)))
            ++ntested;
        }
    REQUIRE(ntested > 0);
  }
#endif
}

// re-initializing an Engine for a task whose AM limit is below the Engine's
// lmax must throw Engine::lmax_exceeded and leave the Engine without kernels:
// the next compute() must throw Engine::missing_kernel rather than use the
// kernel table of the previous task
TEST_CASE("lmax_exceeded leaves no kernels", "[engine][2-body]") {
#if defined(LIBINT2_SUPPORT_ERI) && defined(LIBINT2_SUPPORT_ERI3) && \
    LIBINT2_MAX_AM_3eri > LIBINT2_MAX_AM_eri
  const int lmax = LIBINT2_MAX_AM_3eri;  // valid for 3eri, too high for eri
  const Shell fit{{1.5}, {{lmax, true, {1.0}}}, {{-1.0, 0.5, 0.0}}};
  const Shell s{{1.0}, {{0, false, {1.0}}}, {{0.0, 0.0, 0.0}}};
  const auto eps = std::numeric_limits<double>::epsilon();
  Engine engine(Operator::coulomb, 1, lmax, 0, eps,
                libint2::default_params(Operator::coulomb), BraKet::xs_xx);
  REQUIRE(engine.compute(fit, s, s)[0] != nullptr);
  REQUIRE_THROWS_AS(engine.set(BraKet::xx_xx), Engine::lmax_exceeded);
  CHECK_THROWS_AS(engine.compute(s, s, s, s), Engine::missing_kernel);
  CHECK_THROWS_AS(engine.compute(fit, s, s), Engine::missing_kernel);
#endif
}

TEST_CASE("RKB Coulomb integrals", "[engine][2-body]") {
  const std::vector<Shell> all_shells{
      // pseudorandom s
      Shell{{1.0}, {{0, false, {1.0}}}, {{0.0, 0.0, 0.0}}},
      // pseudorandom p
      Shell{{2.0}, {{1, false, {1.0}}}, {{1.0, 1.0, 1.0}}},
      // pseudorandom contracted p
      Shell{{3.0, 0.7}, {{1, false, {0.4, 0.8}}}, {{-0.5, 0.2, 0.9}}},
      // pseudorandom d
      Shell{{1.2}, {{2, false, {1.0}}}, {{0.3, -0.8, 0.4}}}};
  // the subset of all_shells supported by a kernel with max AM lmax
  auto shells_upto = [&all_shells](int lmax) {
    std::vector<Shell> result;
    for (const auto &sh : all_shells)
      if (sh.contr[0].l <= lmax) result.push_back(sh);
    return result;
  };

  typedef std::array<unsigned int, 12> der_idx;

  // e.g. d_xx maps the derivative index of derivative w.r.t x
  // coord of ket1 and x coord of ket2 in Chemist notation.
  // deriv indices for (LL|SS)
  der_idx d_xx = {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0};
  der_idx d_yy = {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0};
  der_idx d_zz = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1};
  der_idx d_yz = {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1};
  der_idx d_zy = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0};
  der_idx d_zx = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0};
  der_idx d_xz = {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  der_idx d_xy = {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0};
  der_idx d_yx = {0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0};

#if LIBINT2_TASK_EXISTS_4coulomb_opop && LIBINT2_TASK_EXISTS_4opop_coulomb_opop
  SECTION("Coulombσpσp and σpσpCoulombσpσp") {
    // the reference is expensive: skip the uncontracted p shell (the
    // contracted p shell covers p)
    auto obs = shells_upto(std::min(LIBINT2_MAX_AM_coulomb_opop,
                                    LIBINT2_MAX_AM_opop_coulomb_opop));
    obs.erase(std::remove_if(obs.begin(), obs.end(),
                             [](const Shell &sh) {
                               return sh.contr[0].l == 1 && sh.nprim() == 1;
                             }),
              obs.end());
    const auto max_nprim = libint2::max_nprim(obs);
    const auto max_l = libint2::max_l(obs);
    Engine engine_llss(Operator::coulomb_opop, max_nprim, max_l, 0);
    Engine engine_ssss(Operator::opop_coulomb_opop, max_nprim, max_l, 0);

    const auto nshell = obs.size();
    for (int s0 = 0; s0 != nshell; ++s0) {
      for (int s1 = 0; s1 != nshell; ++s1) {
        for (int s2 = 0; s2 != nshell; ++s2) {
          for (int s3 = 0; s3 != nshell; ++s3) {
            // the reference is expensive, include at most 1 d shell
            if ((obs[s0].contr[0].l > 1) + (obs[s1].contr[0].l > 1) +
                    (obs[s2].contr[0].l > 1) + (obs[s3].contr[0].l > 1) >
                1)
              continue;
            const auto &results_llss =
                engine_llss.compute(obs[s0], obs[s1], obs[s2], obs[s3]);
            const auto &results_ssss =
                engine_ssss.compute(obs[s0], obs[s1], obs[s2], obs[s3]);
            // 4 buffers for single-spin quaternion components
            REQUIRE(results_llss.size() == 4);
            // 16 buffers for two independent spin-space quaternions
            REQUIRE(results_ssss.size() == 16);

            LIBINT2_REF_REALTYPE Aref[3];
            for (int i = 0; i < 3; ++i) Aref[i] = obs[s0].O[i];
            LIBINT2_REF_REALTYPE Bref[3];
            for (int i = 0; i < 3; ++i) Bref[i] = obs[s1].O[i];
            LIBINT2_REF_REALTYPE Cref[3];
            for (int i = 0; i < 3; ++i) Cref[i] = obs[s2].O[i];
            LIBINT2_REF_REALTYPE Dref[3];
            for (int i = 0; i < 3; ++i) Dref[i] = obs[s3].O[i];

            int ijkl = 0;

            int l0, m0, n0;
            FOR_CART(l0, m0, n0, obs[s0].contr[0].l)

            int l1, m1, n1;
            FOR_CART(l1, m1, n1, obs[s1].contr[0].l)

            int l2, m2, n2;
            FOR_CART(l2, m2, n2, obs[s2].contr[0].l)

            int l3, m3, n3;
            FOR_CART(l3, m3, n3, obs[s3].contr[0].l)

            std::array<LIBINT2_REF_REALTYPE, 4> ref_coulomb_opop{0.0, 0.0, 0.0,
                                                                 0.0};
            std::array<LIBINT2_REF_REALTYPE, 16> ref_opop_coulomb_opop{};
            ref_opop_coulomb_opop.fill(0.0);
            uint p0123 = 0;
            for (uint p0 = 0; p0 < obs[s0].nprim(); p0++) {
              for (uint p1 = 0; p1 < obs[s1].nprim(); p1++) {
                for (uint p2 = 0; p2 < obs[s2].nprim(); p2++) {
                  for (uint p3 = 0; p3 < obs[s3].nprim(); p3++, p0123++) {
                    const LIBINT2_REF_REALTYPE alpha0 = obs[s0].alpha[p0];
                    const LIBINT2_REF_REALTYPE alpha1 = obs[s1].alpha[p1];
                    const LIBINT2_REF_REALTYPE alpha2 = obs[s2].alpha[p2];
                    const LIBINT2_REF_REALTYPE alpha3 = obs[s3].alpha[p3];

                    const LIBINT2_REF_REALTYPE c0 = obs[s0].contr[0].coeff[p0];
                    const LIBINT2_REF_REALTYPE c1 = obs[s1].contr[0].coeff[p1];
                    const LIBINT2_REF_REALTYPE c2 = obs[s2].contr[0].coeff[p2];
                    const LIBINT2_REF_REALTYPE c3 = obs[s3].contr[0].coeff[p3];
                    const LIBINT2_REF_REALTYPE c0123 = c0 * c1 * c2 * c3;

                    auto eri_drrrr = [&](der_idx d_rrrr) {
                      return eri(d_rrrr.data(), l0, m0, n0, alpha0, Aref, l1,
                                 m1, n1, alpha1, Bref, l2, m2, n2, alpha2, Cref,
                                 l3, m3, n3, alpha3, Dref, 0);
                    };

                    // helper: build der_idx from 4 derivative directions
                    // (0=x, 1=y, 2=z) for centers A, B, C, D
                    constexpr int X = 0, Y = 1, Z = 2;
                    auto didx = [](int a, int b, int c, int d) -> der_idx {
                      der_idx r = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
                      r[a] = 1;
                      r[3 + b] = 1;
                      r[6 + c] = 1;
                      r[9 + d] = 1;
                      return r;
                    };
                    // shorthand: evaluate derivative ERI from 4 directions
                    auto D = [&](int a, int b, int c, int d) {
                      return eri_drrrr(didx(a, b, c, d));
                    };

                    // (LL|SS)
                    ref_coulomb_opop[0] +=
                        c0123 *
                        (eri_drrrr(d_xx) + eri_drrrr(d_yy) + eri_drrrr(d_zz));
                    ref_coulomb_opop[1] +=
                        c0123 * (eri_drrrr(d_yz) - eri_drrrr(d_zy));
                    ref_coulomb_opop[2] +=
                        c0123 * (eri_drrrr(d_zx) - eri_drrrr(d_xz));
                    ref_coulomb_opop[3] +=
                        c0123 * (eri_drrrr(d_xy) - eri_drrrr(d_yx));

                    // (SS|SS) — 16 components, Option A: index = 4*bra + ket
                    // 0:SS  1:SX  2:SY  3:SZ
                    ref_opop_coulomb_opop[0] +=
                        c0123 * (D(X, X, X, X) + D(X, X, Y, Y) + D(X, X, Z, Z) +
                                 D(Y, Y, X, X) + D(Y, Y, Y, Y) + D(Y, Y, Z, Z) +
                                 D(Z, Z, X, X) + D(Z, Z, Y, Y) + D(Z, Z, Z, Z));
                    ref_opop_coulomb_opop[1] +=
                        c0123 * (D(X, X, Y, Z) - D(X, X, Z, Y) + D(Y, Y, Y, Z) -
                                 D(Y, Y, Z, Y) + D(Z, Z, Y, Z) - D(Z, Z, Z, Y));
                    ref_opop_coulomb_opop[2] +=
                        c0123 * (D(X, X, Z, X) - D(X, X, X, Z) + D(Y, Y, Z, X) -
                                 D(Y, Y, X, Z) + D(Z, Z, Z, X) - D(Z, Z, X, Z));
                    ref_opop_coulomb_opop[3] +=
                        c0123 * (D(X, X, X, Y) - D(X, X, Y, X) + D(Y, Y, X, Y) -
                                 D(Y, Y, Y, X) + D(Z, Z, X, Y) - D(Z, Z, Y, X));
                    // 4:XS  5:XX  6:XY  7:XZ
                    ref_opop_coulomb_opop[4] +=
                        c0123 * (D(Y, Z, X, X) - D(Z, Y, X, X) + D(Y, Z, Y, Y) -
                                 D(Z, Y, Y, Y) + D(Y, Z, Z, Z) - D(Z, Y, Z, Z));
                    ref_opop_coulomb_opop[5] +=
                        c0123 * (-D(Y, Z, Y, Z) + D(Y, Z, Z, Y) +
                                 D(Z, Y, Y, Z) - D(Z, Y, Z, Y));
                    ref_opop_coulomb_opop[6] +=
                        c0123 * (-D(Y, Z, Z, X) + D(Y, Z, X, Z) +
                                 D(Z, Y, Z, X) - D(Z, Y, X, Z));
                    ref_opop_coulomb_opop[7] +=
                        c0123 * (-D(Y, Z, X, Y) + D(Y, Z, Y, X) +
                                 D(Z, Y, X, Y) - D(Z, Y, Y, X));
                    // 8:YS  9:YX  10:YY  11:YZ
                    ref_opop_coulomb_opop[8] +=
                        c0123 * (D(Z, X, X, X) - D(X, Z, X, X) + D(Z, X, Y, Y) -
                                 D(X, Z, Y, Y) + D(Z, X, Z, Z) - D(X, Z, Z, Z));
                    ref_opop_coulomb_opop[9] +=
                        c0123 * (-D(Z, X, Y, Z) + D(Z, X, Z, Y) +
                                 D(X, Z, Y, Z) - D(X, Z, Z, Y));
                    ref_opop_coulomb_opop[10] +=
                        c0123 * (-D(Z, X, Z, X) + D(Z, X, X, Z) +
                                 D(X, Z, Z, X) - D(X, Z, X, Z));
                    ref_opop_coulomb_opop[11] +=
                        c0123 * (-D(Z, X, X, Y) + D(Z, X, Y, X) +
                                 D(X, Z, X, Y) - D(X, Z, Y, X));
                    // 12:ZS  13:ZX  14:ZY  15:ZZ
                    ref_opop_coulomb_opop[12] +=
                        c0123 * (D(X, Y, X, X) - D(Y, X, X, X) + D(X, Y, Y, Y) -
                                 D(Y, X, Y, Y) + D(X, Y, Z, Z) - D(Y, X, Z, Z));
                    ref_opop_coulomb_opop[13] +=
                        c0123 * (-D(X, Y, Y, Z) + D(X, Y, Z, Y) +
                                 D(Y, X, Y, Z) - D(Y, X, Z, Y));
                    ref_opop_coulomb_opop[14] +=
                        c0123 * (-D(X, Y, Z, X) + D(X, Y, X, Z) +
                                 D(Y, X, Z, X) - D(Y, X, X, Z));
                    ref_opop_coulomb_opop[15] +=
                        c0123 * (-D(X, Y, X, Y) + D(X, Y, Y, X) +
                                 D(Y, X, X, Y) - D(Y, X, Y, X));
                  }
                }
              }
            }

            const double ABSOLUTE_DEVIATION_THRESHOLD = 5.0E-14;
            const double RELATIVE_DEVIATION_THRESHOLD =
                1.0E-9;  // For more detail on choice of these thresholds, see
                         // the comments in the TEST_CASE "eri geometric
                         // derivatives"

            std::array<LIBINT2_REF_REALTYPE, 4> abs_errs_llss;
            std::array<LIBINT2_REF_REALTYPE, 4> rel_abs_errs_llss;

            // (LL|SS) has 4 components
            for (auto comp = 0; comp < 4; ++comp) {
              abs_errs_llss[comp] =
                  abs(ref_coulomb_opop[comp] - results_llss[comp][ijkl]);
              rel_abs_errs_llss[comp] =
                  abs(abs_errs_llss[comp] / ref_coulomb_opop[comp]);

              bool llss_not_ok =
                  rel_abs_errs_llss[comp] > RELATIVE_DEVIATION_THRESHOLD &&
                  abs_errs_llss[comp] > ABSOLUTE_DEVIATION_THRESHOLD;

              if (llss_not_ok) {
                std::cout << "(l0 l1| l2 l3) = "
                          << "(" << s0 << " " << s1 << " | " << s2 << " " << s3
                          << ") "
                          << "Elem " << ijkl << " comp= " << comp
                          << " : ref = " << ref_coulomb_opop[comp]
                          << " libint = " << results_llss[comp][ijkl]
                          << " relabs_error = " << rel_abs_errs_llss[comp]
                          << " abs_error = " << abs_errs_llss[comp]
                          << std::endl;
              }
              REQUIRE(!llss_not_ok);
            }

            // (SS|SS) has 16 components (two independent spin spaces); the
            // reference involves 4th derivatives and loses more precision
            const double SSSS_ABSOLUTE_DEVIATION_THRESHOLD = 1.0E-12;
            for (auto comp = 0; comp < 16; ++comp) {
              auto abs_err_ssss =
                  abs(ref_opop_coulomb_opop[comp] - results_ssss[comp][ijkl]);
              auto rel_abs_err_ssss =
                  abs(abs_err_ssss / ref_opop_coulomb_opop[comp]);

              bool ssss_not_ok =
                  rel_abs_err_ssss > RELATIVE_DEVIATION_THRESHOLD &&
                  abs_err_ssss > SSSS_ABSOLUTE_DEVIATION_THRESHOLD;

              if (ssss_not_ok) {
                std::cout << "(l0 l1| l2 l3) = "
                          << "(" << s0 << " " << s1 << " | " << s2 << " " << s3
                          << ") "
                          << "Elem " << ijkl << " comp= " << comp
                          << " : ref = " << ref_opop_coulomb_opop[comp]
                          << " libint = " << results_ssss[comp][ijkl]
                          << " relabs_error = " << rel_abs_err_ssss
                          << " abs_error = " << abs_err_ssss << std::endl;
              }
              REQUIRE(!ssss_not_ok);
            }

            ++ijkl;
            END_FOR_CART
            END_FOR_CART
            END_FOR_CART
            END_FOR_CART
          }
        }
      }
    }
  }

#endif  // 4-center Coulombσpσp and σpσpCoulombσpσp

#if LIBINT2_TASK_EXISTS_4op_coulomb_op
  SECTION("op_coulomb_op") {
    const auto obs = shells_upto(LIBINT2_MAX_AM_op_coulomb_op);
    const auto max_nprim = libint2::max_nprim(obs);
    const auto max_l = libint2::max_l(obs);
    Engine engine_opCop(Operator::op_coulomb_op, max_nprim, max_l, 0);

    const auto nshell = obs.size();
    for (int s0 = 0; s0 != nshell; ++s0) {
      for (int s1 = 0; s1 != nshell; ++s1) {
        for (int s2 = 0; s2 != nshell; ++s2) {
          for (int s3 = 0; s3 != nshell; ++s3) {
            const auto &results =
                engine_opCop.compute(obs[s0], obs[s1], obs[s2], obs[s3]);
            REQUIRE(results.size() == 9);

            LIBINT2_REF_REALTYPE Aref[3], Bref[3], Cref[3], Dref[3];
            for (int i = 0; i < 3; ++i) Aref[i] = obs[s0].O[i];
            for (int i = 0; i < 3; ++i) Bref[i] = obs[s1].O[i];
            for (int i = 0; i < 3; ++i) Cref[i] = obs[s2].O[i];
            for (int i = 0; i < 3; ++i) Dref[i] = obs[s3].O[i];

            int ijkl = 0;

            int l0, m0, n0;
            FOR_CART(l0, m0, n0, obs[s0].contr[0].l)
            int l1, m1, n1;
            FOR_CART(l1, m1, n1, obs[s1].contr[0].l)
            int l2, m2, n2;
            FOR_CART(l2, m2, n2, obs[s2].contr[0].l)
            int l3, m3, n3;
            FOR_CART(l3, m3, n3, obs[s3].contr[0].l)

            // Raw 3x3 dyadic T_{ab} = (a b∂_a | c d∂_b), accumulated in
            // chemist-notation index 3*a+b. After the primitive loop we
            // project into the 9 SO(3) irreducible components that the engine
            // returns: Scalar, AntisymX/Y/Z, SymTLDiagA/B, SymTLOffXY/XZ/YZ.
            std::array<LIBINT2_REF_REALTYPE, 9> ref_raw{};
            ref_raw.fill(0.0);

            for (uint p0 = 0; p0 < obs[s0].nprim(); p0++) {
              for (uint p1 = 0; p1 < obs[s1].nprim(); p1++) {
                for (uint p2 = 0; p2 < obs[s2].nprim(); p2++) {
                  for (uint p3 = 0; p3 < obs[s3].nprim(); p3++) {
                    const LIBINT2_REF_REALTYPE alpha0 = obs[s0].alpha[p0];
                    const LIBINT2_REF_REALTYPE alpha1 = obs[s1].alpha[p1];
                    const LIBINT2_REF_REALTYPE alpha2 = obs[s2].alpha[p2];
                    const LIBINT2_REF_REALTYPE alpha3 = obs[s3].alpha[p3];
                    const LIBINT2_REF_REALTYPE c0 = obs[s0].contr[0].coeff[p0];
                    const LIBINT2_REF_REALTYPE c1 = obs[s1].contr[0].coeff[p1];
                    const LIBINT2_REF_REALTYPE c2 = obs[s2].contr[0].coeff[p2];
                    const LIBINT2_REF_REALTYPE c3 = obs[s3].contr[0].coeff[p3];
                    const LIBINT2_REF_REALTYPE c0123 = c0 * c1 * c2 * c3;

                    // Deriv on ν (center B, index 1) and on λ (center D, idx
                    // 3).
                    auto didx_bd = [](int a, int b) -> der_idx {
                      der_idx r = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
                      r[3 + a] = 1;
                      r[9 + b] = 1;
                      return r;
                    };
                    auto D = [&](int a, int b) {
                      auto di = didx_bd(a, b);
                      return eri(di.data(), l0, m0, n0, alpha0, Aref, l1, m1,
                                 n1, alpha1, Bref, l2, m2, n2, alpha2, Cref, l3,
                                 m3, n3, alpha3, Dref, 0);
                    };
                    for (int a = 0; a < 3; ++a)
                      for (int b = 0; b < 3; ++b)
                        ref_raw[3 * a + b] += c0123 * D(a, b);
                  }
                }
              }
            }

            // Project raw dyadic into the 9 SO(3) irrep components used by
            // op_coulomb_op (must match σpCoulombσp_Descr::Component order).
            const auto Txx = ref_raw[0];
            const auto Txy = ref_raw[1];
            const auto Txz = ref_raw[2];
            const auto Tyx = ref_raw[3];
            const auto Tyy = ref_raw[4];
            const auto Tyz = ref_raw[5];
            const auto Tzx = ref_raw[6];
            const auto Tzy = ref_raw[7];
            const auto Tzz = ref_raw[8];
            std::array<LIBINT2_REF_REALTYPE, 9> ref_op_coulomb_op{
                Txx + Tyy + Tzz,        // Scalar
                Tyz - Tzy,              // AntisymX
                Tzx - Txz,              // AntisymY
                Txy - Tyx,              // AntisymZ
                Txx - Tyy,              // SymTLDiagA
                2.0 * Tzz - Txx - Tyy,  // SymTLDiagB
                Txy + Tyx,              // SymTLOffXY
                Txz + Tzx,              // SymTLOffXZ
                Tyz + Tzy,              // SymTLOffYZ
            };

            const double ABSOLUTE_DEVIATION_THRESHOLD = 5.0E-14;
            const double RELATIVE_DEVIATION_THRESHOLD = 1.0E-9;
            for (auto comp = 0; comp < 9; ++comp) {
              auto abs_err = abs(ref_op_coulomb_op[comp] - results[comp][ijkl]);
              auto rel_abs_err = abs(abs_err / ref_op_coulomb_op[comp]);
              bool not_ok = rel_abs_err > RELATIVE_DEVIATION_THRESHOLD &&
                            abs_err > ABSOLUTE_DEVIATION_THRESHOLD;
              if (not_ok) {
                std::cout << "(l0 l1| l2 l3) = (" << s0 << " " << s1 << " | "
                          << s2 << " " << s3 << ") Elem " << ijkl
                          << " comp= " << comp
                          << " : ref = " << ref_op_coulomb_op[comp]
                          << " libint = " << results[comp][ijkl]
                          << " relabs_error = " << rel_abs_err
                          << " abs_error = " << abs_err << std::endl;
              }
              REQUIRE(!not_ok);
            }

            ++ijkl;
            END_FOR_CART
            END_FOR_CART
            END_FOR_CART
            END_FOR_CART
          }
        }
      }
    }
  }

#endif  // 4-center op_coulomb_op

// The 3-center RKB sections below exercise kernels that exist only when
// LIBINT2_ENABLE_RKB_ERI3 >= 0 (else compute() throws Engine::missing_kernel).
// Guard at compile time on the generated task-exists macros.
#if LIBINT2_TASK_EXISTS_3coulomb_opop && LIBINT2_TASK_EXISTS_4coulomb_opop
  SECTION("Coulombσpσp 3-center xs_xx") {
    // 3-center RKB Coulombσpσp integral, (P | σ·p_μ σ·p_ν / r12) with P on
    // the bra (fitting/DF) and σ·p acting on the AO pair (μ,ν) on the ket.
    //
    // Reference: the same operator and BraKet computed via the 4-center
    // engine with Shell::unit() at the dummy bra position. The 4-center
    // engine path is already validated by the "Coulombσpσp and
    // σpσpCoulombσpσp" SECTION above against a primitive-eri reference, so
    // by transitivity the 4-center engine is a trustworthy reference for the
    // 3-center engine. Both code paths share the dummy-shell trick
    // internally, but they are independently generated (different task
    // labels, different dispatch tables), so this catches dispatch-table /
    // codegen wiring mistakes in the new 3-center path.

    Shell dfsh{{1.5}, {{0, false, {1.0}}}, {{-1.0, 0.5, 0.0}}};
    std::vector<Shell> dfs{dfsh};
    const auto &unitshell = libint2::Shell::unit();

    const auto obs = shells_upto(
        std::min(LIBINT2_MAX_AM_3coulomb_opop, LIBINT2_MAX_AM_coulomb_opop));
    const auto max_nprim = libint2::max_nprim(obs);
    const auto max_l = libint2::max_l(obs);
    Engine engine_3c(Operator::coulomb_opop,
                     std::max(max_nprim, libint2::max_nprim(dfs)),
                     std::max(max_l, libint2::max_l(dfs)), 0);
    engine_3c.set(BraKet::xs_xx);
    // 4-center reference — same operator, default BraKet::xx_xx
    Engine engine_4c_ref(Operator::coulomb_opop,
                         std::max(max_nprim, libint2::max_nprim(dfs)),
                         std::max(max_l, libint2::max_l(dfs)), 0);

    const auto nshell = obs.size();
    for (int sa = 0; sa != nshell; ++sa) {
      for (int sb = 0; sb != nshell; ++sb) {
        const auto &results_3c = engine_3c.compute(dfsh, obs[sa], obs[sb]);
        // 4-center reference: (dfsh, unit | sh_a, sh_b). σ·p still acts on
        // (C, D) = (sh_a, sh_b), matching the 3-center xs_xx mapping.
        const auto &results_4c_ref =
            engine_4c_ref.compute(dfsh, unitshell, obs[sa], obs[sb]);
        REQUIRE(results_3c.size() == 4);
        REQUIRE(results_4c_ref.size() == 4);

        const auto n_df = dfsh.size();
        const auto n_a = obs[sa].size();
        const auto n_b = obs[sb].size();
        const auto n_total = n_df * n_a * n_b;
        // 4-center buffer is (df, unit, sa, sb); n_unit = 1; so the linear
        // layout coincides with the 3-center (df, sa, sb) layout.

        const double ABS_TOL = 5.0E-14;
        const double REL_TOL = 1.0E-9;
        for (auto comp = 0; comp < 4; ++comp) {
          for (size_t k = 0; k < n_total; ++k) {
            const auto v_3c = results_3c[comp][k];
            const auto v_ref = results_4c_ref[comp][k];
            const auto abs_err = std::abs(v_3c - v_ref);
            const auto rel_err =
                std::abs(v_ref) > 1e-30 ? std::abs(abs_err / v_ref) : 0.0;
            bool not_ok = rel_err > REL_TOL && abs_err > ABS_TOL;
            if (not_ok) {
              std::cout << "(df | sa=" << sa << " sb=" << sb
                        << ") comp=" << comp << " elem=" << k << ": 3c=" << v_3c
                        << " 4c_ref=" << v_ref << " abs_err=" << abs_err
                        << " rel_err=" << rel_err << std::endl;
            }
            REQUIRE(!not_ok);
          }
        }
      }
    }
  }

  SECTION("Coulombσpσp 3-center xx_xs alias") {
    // The xx_xs braket re-routes through the same xs_xx kernel via the
    // Engine::compute2 bra↔ket swap. Verify the alias yields the same
    // integral values, with the output buffer in the user-requested
    // (sh_a, sh_b | dfsh) index layout.

    Shell dfsh{{1.5}, {{0, false, {1.0}}}, {{-1.0, 0.5, 0.0}}};
    std::vector<Shell> dfs{dfsh};

    const auto obs = shells_upto(LIBINT2_MAX_AM_3coulomb_opop);
    const auto max_nprim = libint2::max_nprim(obs);
    const auto max_l = libint2::max_l(obs);
    Engine eng_xs(Operator::coulomb_opop,
                  std::max(max_nprim, libint2::max_nprim(dfs)),
                  std::max(max_l, libint2::max_l(dfs)), 0);
    eng_xs.set(BraKet::xs_xx);
    Engine eng_xxs(Operator::coulomb_opop,
                   std::max(max_nprim, libint2::max_nprim(dfs)),
                   std::max(max_l, libint2::max_l(dfs)), 0);
    eng_xxs.set(BraKet::xx_xs);

    for (size_t sa = 0; sa != obs.size(); ++sa) {
      for (size_t sb = 0; sb != obs.size(); ++sb) {
        const auto n_df = dfsh.size();
        const auto n_a = obs[sa].size();
        const auto n_b = obs[sb].size();

        const auto &res_xs = eng_xs.compute(dfsh, obs[sa], obs[sb]);
        // re-set xs braket each loop in case Catch2 SECTION re-entry resets
        eng_xxs.set(BraKet::xx_xs);
        const auto &res_xxs = eng_xxs.compute(obs[sa], obs[sb], dfsh);

        for (auto c = 0; c < 4; ++c) {
          for (size_t i_df = 0; i_df < n_df; ++i_df) {
            for (size_t i_a = 0; i_a < n_a; ++i_a) {
              for (size_t i_b = 0; i_b < n_b; ++i_b) {
                const auto v_xs = res_xs[c][i_df * n_a * n_b + i_a * n_b + i_b];
                const auto v_xxs =
                    res_xxs[c][i_a * n_b * n_df + i_b * n_df + i_df];
                const auto abs_err = abs(v_xs - v_xxs);
                const auto rel_err =
                    std::abs(v_xs) > 1e-30 ? double(abs(abs_err / v_xs)) : 0.0;
                bool not_ok = rel_err > 1.0E-9 && abs_err > 5.0E-14;
                if (not_ok) {
                  std::cout
                      << "xx_xs alias mismatch: (sa=" << sa << " sb=" << sb
                      << ") comp=" << c << " idx(df=" << i_df << ",a=" << i_a
                      << ",b=" << i_b << "): xs=" << v_xs << " xxs=" << v_xxs
                      << " abs_err=" << abs_err << std::endl;
                }
                REQUIRE(!not_ok);
              }
            }
          }
        }
      }
    }
  }
#endif  // 3-center Coulombσpσp

#if LIBINT2_TASK_EXISTS_3coulomb_op
  SECTION("Coulombσp 3-center xs_xx") {
    // 3-center single-σ·p integral
    //   (P | 1/r12 | μ, σ·p ν),
    // with the fitting function P on the bra (BraKet::xs_xx) and a single σ·p
    // acting on the 2nd ket AO function ν only. The 3 Cartesian components
    // (x,y,z) are the bare gradient of ν, so the reference is the 3-center ERI
    // with a *single first derivative on ν* (the shell that carries σ·p) — the
    // other geometric derivatives are not needed.
    //
    // The 3-center ERI is evaluated with the primitive eri() reference by
    // placing a unit s-shell (alpha=0, at the origin) at the dummy bra-2
    // position, exactly as the engine does internally via Shell::unit(): the
    // 3-arg compute(df, μ, ν) maps to the 4-center order (df, unit, μ, ν) =
    // (A,B,C,D), so σ·p on ν is the derivative on center D (der_idx slots
    // 9..11). See the "Coulombσpσp and σpσpCoulombσpσp" SECTION for the
    // (un-derivatived) primitive-eri reference pattern.

    Shell dfsh{{1.5}, {{0, false, {1.0}}}, {{-1.0, 0.5, 0.0}}};
    std::vector<Shell> dfs{dfsh};

    const auto obs = shells_upto(LIBINT2_MAX_AM_3coulomb_op);
    const auto max_nprim = libint2::max_nprim(obs);
    const auto max_l = libint2::max_l(obs);
    Engine engine_3c(Operator::coulomb_op,
                     std::max(max_nprim, libint2::max_nprim(dfs)),
                     std::max(max_l, libint2::max_l(dfs)), 0);
    engine_3c.set(BraKet::xs_xx);
    // coulomb_op has no 4-center kernel, hence default_braket(coulomb_op)
    // must be the 3-center xs_xx, so engine_3c_default — built WITHOUT an
    // explicit set() — must produce results identical to engine_3c.
    CHECK(default_braket(Operator::coulomb_op) == BraKet::xs_xx);
    Engine engine_3c_default(Operator::coulomb_op,
                             std::max(max_nprim, libint2::max_nprim(dfs)),
                             std::max(max_l, libint2::max_l(dfs)), 0);
    // no 4-center coulomb_op kernel -> must throw rather than silently
    // computing plain Coulomb integrals
    {
      Engine engine_4c(Operator::coulomb_op,
                       std::max(max_nprim, libint2::max_nprim(dfs)),
                       std::max(max_l, libint2::max_l(dfs)), 0);
      engine_4c.set(BraKet::xx_xx);
      CHECK_THROWS_AS(engine_4c.compute(obs[0], obs[0], obs[0], obs[0]),
                      Engine::missing_kernel);
    }

    // unit (dummy) shell at the bra-2 position: alpha=0, l=0, at the origin.
    LIBINT2_REF_REALTYPE Bref[3] = {0.0, 0.0, 0.0};
    const LIBINT2_REF_REALTYPE alpha_unit = 0.0;

    const auto nshell = obs.size();
    for (int sa = 0; sa != nshell; ++sa) {
      for (int sb = 0; sb != nshell; ++sb) {
        const auto &results_3c = engine_3c.compute(dfsh, obs[sa], obs[sb]);
        REQUIRE(results_3c.size() == 3);
        // default-braket engine (no set()) must match the explicit-xs_xx engine
        const auto &results_3c_default =
            engine_3c_default.compute(dfsh, obs[sa], obs[sb]);
        REQUIRE(results_3c_default.size() == 3);

        LIBINT2_REF_REALTYPE Aref[3], Cref[3], Dref[3];
        for (int i = 0; i < 3; ++i) Aref[i] = dfsh.O[i];
        for (int i = 0; i < 3; ++i) Cref[i] = obs[sa].O[i];
        for (int i = 0; i < 3; ++i) Dref[i] = obs[sb].O[i];

        int ijk = 0;

        int l0, m0, n0;
        FOR_CART(l0, m0, n0, dfsh.contr[0].l)  // P (fitting, s)
        int l2, m2, n2;
        FOR_CART(l2, m2, n2, obs[sa].contr[0].l)  // μ
        int l3, m3, n3;
        FOR_CART(l3, m3, n3, obs[sb].contr[0].l)  // ν (σ·p acts here)

        std::array<LIBINT2_REF_REALTYPE, 3> ref_coulomb_op{0.0, 0.0, 0.0};

        for (uint p0 = 0; p0 < dfsh.nprim(); p0++) {
          for (uint p2 = 0; p2 < obs[sa].nprim(); p2++) {
            for (uint p3 = 0; p3 < obs[sb].nprim(); p3++) {
              const LIBINT2_REF_REALTYPE alpha0 = dfsh.alpha[p0];
              const LIBINT2_REF_REALTYPE alpha2 = obs[sa].alpha[p2];
              const LIBINT2_REF_REALTYPE alpha3 = obs[sb].alpha[p3];
              const LIBINT2_REF_REALTYPE c0 = dfsh.contr[0].coeff[p0];
              const LIBINT2_REF_REALTYPE c2 = obs[sa].contr[0].coeff[p2];
              const LIBINT2_REF_REALTYPE c3 = obs[sb].contr[0].coeff[p3];
              const LIBINT2_REF_REALTYPE c023 = c0 * c2 * c3;  // unit coeff = 1

              // single first derivative on ν (center D), direction k.
              auto eri_dD = [&](int k) {
                der_idx d = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
                d[9 + k] = 1;  // slots 9,10,11 = D_x, D_y, D_z
                return eri(d.data(), l0, m0, n0, alpha0, Aref, 0, 0, 0,
                           alpha_unit, Bref, l2, m2, n2, alpha2, Cref, l3, m3,
                           n3, alpha3, Dref, 0);
              };

              for (int k = 0; k < 3; ++k) ref_coulomb_op[k] += c023 * eri_dD(k);
            }
          }
        }

        const double ABS_TOL = 5.0E-14;
        const double REL_TOL = 1.0E-9;
        for (auto comp = 0; comp < 3; ++comp) {
          const auto v_3c = results_3c[comp][ijk];
          CHECK(results_3c_default[comp][ijk] == v_3c);
          const auto v_ref = ref_coulomb_op[comp];
          const auto abs_err = std::abs(v_3c - v_ref);
          const auto rel_err =
              std::abs(v_ref) > 1e-30 ? std::abs(abs_err / v_ref) : 0.0;
          bool not_ok = rel_err > REL_TOL && abs_err > ABS_TOL;
          if (not_ok) {
            std::cout << "(df | sa=" << sa << " sb=" << sb << ") comp=" << comp
                      << " elem=" << ijk << ": 3c=" << v_3c << " ref=" << v_ref
                      << " abs_err=" << abs_err << " rel_err=" << rel_err
                      << std::endl;
          }
          REQUIRE(!not_ok);
        }

        ++ijk;
        END_FOR_CART
        END_FOR_CART
        END_FOR_CART
      }
    }
  }
#endif  // LIBINT2_TASK_EXISTS_3coulomb_op
}

// RKB integrals over solid harmonics shells. The RKB reference tests above use
// cartesian shells only; the engine's cartesian-to-solid-harmonics transform
// shares the scratch/permute path with the per-component sign corrections of
// the RKB operators, so verify here that the two compose: the result for pure
// shells must equal the cartesian result transformed with the library's own
// solid harmonics routines. Also covers the xx_xs alias of coulomb_op and a
// solid harmonics fitting center.
namespace {

typedef std::vector<std::vector<double>> rkb_shellsets_t;

rkb_shellsets_t rkb_copy_results(const Engine::target_ptr_vec &results,
                                 const std::vector<Shell> &shells) {
  size_t n = 1;
  for (const auto &sh : shells) n *= sh.size();
  rkb_shellsets_t sets;
  for (const auto *target : results) {
    REQUIRE(target != nullptr);
    sets.emplace_back(target, target + n);
  }
  return sets;
}

/// transforms shell sets computed over @c src_shells to the basis of
/// @c tgt_shells; the two differ only in which shells are pure
rkb_shellsets_t rkb_to_pure(const rkb_shellsets_t &sets,
                            const std::vector<Shell> &src_shells,
                            const std::vector<Shell> &tgt_shells) {
  REQUIRE(src_shells.size() == tgt_shells.size());
  const auto nsh = src_shells.size();
  rkb_shellsets_t result;
  for (const auto &set : sets) {
    auto src = set;
    for (size_t k = 0; k != nsh; ++k) {
      const auto l = src_shells[k].contr[0].l;
      REQUIRE(l == tgt_shells[k].contr[0].l);
      if (src_shells[k].contr[0].pure == tgt_shells[k].contr[0].pure || l == 0)
        continue;
      REQUIRE(tgt_shells[k].contr[0].pure);
      // dimensions before k are already transformed, those after are not
      size_t n1 = 1, n2 = 1;
      for (size_t i = 0; i != k; ++i) n1 *= tgt_shells[i].size();
      for (size_t i = k + 1; i != nsh; ++i) n2 *= src_shells[i].size();
      std::vector<double> tgt(n1 * (2 * l + 1) * n2);
      if (k == 0)
        libint2::solidharmonics::transform_first(l, n2, src.data(), tgt.data());
      else if (k + 1 == nsh)
        libint2::solidharmonics::transform_last(n1, l, src.data(), tgt.data());
      else
        libint2::solidharmonics::transform_inner(n1, l, n2, src.data(),
                                                 tgt.data());
      src = std::move(tgt);
    }
    result.push_back(std::move(src));
  }
  return result;
}

void rkb_require_equal(const rkb_shellsets_t &a, const rkb_shellsets_t &b,
                       const std::string &what) {
  REQUIRE(a.size() == b.size());
  for (size_t c = 0; c != a.size(); ++c) {
    REQUIRE(a[c].size() == b[c].size());
    for (size_t i = 0; i != a[c].size(); ++i) {
      const auto abs_err = std::abs(a[c][i] - b[c][i]);
      const auto rel_err =
          std::abs(b[c][i]) > 1e-30 ? abs_err / std::abs(b[c][i]) : 0.0;
      INFO(what << " comp=" << c << " elem=" << i << " value=" << a[c][i]
                << " ref=" << b[c][i]);
      REQUIRE(!(rel_err > 1e-10 && abs_err > 1e-13));
    }
  }
}

Shell rkb_shell(const Shell &sh, bool pure) {
  auto result = sh;
  result.contr[0].pure = pure;
  return result;
}

}  // namespace

TEST_CASE("RKB integrals with solid harmonics shells", "[engine][2-body]") {
  // this TU includes test_eri/eri.h, which declares libint2::parity; make sure
  // the solid harmonics coefficients computed from this TU's instantiation
  // (which wins at link time on some toolchains) do not pick it up
  REQUIRE(libint2::solidharmonics::SolidHarmonicsCoefficients<double>::coeff(
              2, 0, 2, 0, 0) == Approx(-0.5));
  // cartesian s, p, contracted p, d; the solid harmonics variants of the same
  const std::vector<Shell> cart_shells{
      Shell{{1.0}, {{0, false, {1.0}}}, {{0.0, 0.0, 0.0}}},
      Shell{{2.0}, {{1, false, {1.0}}}, {{1.0, 1.0, 1.0}}},
      Shell{{3.0, 0.7}, {{1, false, {0.4, 0.8}}}, {{-0.5, 0.2, 0.9}}},
      Shell{{1.2}, {{2, false, {1.0}}}, {{0.3, -0.8, 0.4}}}};
  std::vector<Shell> pure_shells;
  for (const auto &sh : cart_shells) pure_shells.push_back(rkb_shell(sh, true));
  const auto nprim = libint2::max_nprim(cart_shells);

  auto upto = [](const std::vector<Shell> &shells, int lmax) {
    std::vector<Shell> result;
    for (const auto &sh : shells)
      if (sh.contr[0].l <= lmax) result.push_back(sh);
    return result;
  };

  // 4-center: compare pure vs transformed cartesian for every quartet
  auto check_4center = [&](Operator oper, int lmax_oper) {
    const auto cart = upto(cart_shells, lmax_oper);
    const auto pure = upto(pure_shells, lmax_oper);
    Engine engine(oper, nprim, libint2::max_l(cart), 0);
    const auto nsh = cart.size();
    for (size_t a = 0; a != nsh; ++a)
      for (size_t b = 0; b != nsh; ++b)
        for (size_t c = 0; c != nsh; ++c)
          for (size_t d = 0; d != nsh; ++d) {
            const std::vector<Shell> cs{cart[a], cart[b], cart[c], cart[d]};
            const std::vector<Shell> ps{pure[a], pure[b], pure[c], pure[d]};
            const auto ref =
                rkb_to_pure(rkb_copy_results(
                                engine.compute(cs[0], cs[1], cs[2], cs[3]), cs),
                            cs, ps);
            const auto val = rkb_copy_results(
                engine.compute(ps[0], ps[1], ps[2], ps[3]), ps);
            std::ostringstream oss;
            oss << "oper=" << static_cast<int>(oper) << " (" << a << b << "|"
                << c << d << ")";
            rkb_require_equal(val, ref, oss.str());
          }
  };

#if LIBINT2_TASK_EXISTS_4coulomb_opop
  SECTION("coulomb_opop 4-center") {
    check_4center(Operator::coulomb_opop, LIBINT2_MAX_AM_coulomb_opop);
  }
#endif
#if LIBINT2_TASK_EXISTS_4opop_coulomb_opop
  SECTION("opop_coulomb_opop 4-center") {
    check_4center(Operator::opop_coulomb_opop,
                  LIBINT2_MAX_AM_opop_coulomb_opop);
  }
#endif
#if LIBINT2_TASK_EXISTS_4op_coulomb_op
  SECTION("op_coulomb_op 4-center") {
    check_4center(Operator::op_coulomb_op, LIBINT2_MAX_AM_op_coulomb_op);
  }
#endif

  // 3-center: fitting shells s and d, AO shells s, p, d, cartesian vs pure,
  // both brakets
  auto check_3center = [&](Operator oper, int lmax_fit, int lmax_ao) {
    const std::vector<Shell> fit_pure{
        Shell{{1.5}, {{0, true, {1.0}}}, {{-1.0, 0.5, 0.0}}},
        Shell{{0.9}, {{2, true, {1.0}}}, {{-1.0, 0.5, 0.0}}}};
    const auto fit = upto(fit_pure, lmax_fit);
    const auto cart = upto(cart_shells, lmax_ao);
    const auto pure = upto(pure_shells, lmax_ao);
    const auto nprim_3c = std::max(nprim, libint2::max_nprim(fit));
    const auto lmax_3c = std::max(libint2::max_l(cart), libint2::max_l(fit));
    const auto eps = std::numeric_limits<double>::epsilon();
    Engine eng_xs(oper, nprim_3c, lmax_3c, 0, eps,
                  libint2::default_params(oper), BraKet::xs_xx);
    Engine eng_xxs(oper, nprim_3c, lmax_3c, 0, eps,
                   libint2::default_params(oper), BraKet::xx_xs);

    for (size_t p = 0; p != fit.size(); ++p) {
      const auto fit_ref = rkb_shell(fit[p], false);
      for (size_t a = 0; a != cart.size(); ++a)
        for (size_t b = 0; b != cart.size(); ++b) {
          const std::vector<Shell> rs{fit_ref, cart[a], cart[b]};
          const std::vector<Shell> ps{fit[p], pure[a], pure[b]};
          const auto ref = rkb_to_pure(
              rkb_copy_results(eng_xs.compute(rs[0], rs[1], rs[2]), rs), rs,
              ps);
          const auto val =
              rkb_copy_results(eng_xs.compute(ps[0], ps[1], ps[2]), ps);
          std::ostringstream oss;
          oss << "oper=" << static_cast<int>(oper) << " (" << p << "|" << a << b
              << ")";
          rkb_require_equal(val, ref, oss.str() + " xs_xx");

          // xx_xs alias: same integrals, (a b | fit) layout
          const auto val_xxs =
              rkb_copy_results(eng_xxs.compute(ps[1], ps[2], ps[0]),
                               std::vector<Shell>{ps[1], ps[2], ps[0]});
          const auto n_fit = ps[0].size(), n_a = ps[1].size(),
                     n_b = ps[2].size();
          rkb_shellsets_t val_xxs_as_xs(val_xxs.size());
          for (size_t c = 0; c != val_xxs.size(); ++c) {
            val_xxs_as_xs[c].resize(n_fit * n_a * n_b);
            for (size_t i_f = 0; i_f != n_fit; ++i_f)
              for (size_t i_a = 0; i_a != n_a; ++i_a)
                for (size_t i_b = 0; i_b != n_b; ++i_b)
                  val_xxs_as_xs[c][(i_f * n_a + i_a) * n_b + i_b] =
                      val_xxs[c][(i_a * n_b + i_b) * n_fit + i_f];
          }
          rkb_require_equal(val_xxs_as_xs, val, oss.str() + " xx_xs");
        }
    }
  };

#if LIBINT2_TASK_EXISTS_3coulomb_opop
  SECTION("coulomb_opop 3-center") {
    check_3center(Operator::coulomb_opop, LIBINT2_MAX_AM_3coulomb_opop,
                  LIBINT2_MAX_AM_default);
  }
#if LIBINT2_TASK_EXISTS_4coulomb_opop
  SECTION("coulomb_opop 3-center vs 4-center, pure fitting shell") {
    // the 3-center kernel is generated with a dummy shell next to the fitting
    // shell; check against the 4-center kernel, which has no such structure
    const Shell fit{{0.9}, {{2, true, {1.0}}}, {{-1.0, 0.5, 0.0}}};
    const auto lmax_ao = std::min(
        LIBINT2_MAX_AM_default,
        std::min(LIBINT2_MAX_AM_3coulomb_opop, LIBINT2_MAX_AM_coulomb_opop));
    if (LIBINT2_MAX_AM_3coulomb_opop >= 2 && LIBINT2_MAX_AM_coulomb_opop >= 2) {
      const auto pure = upto(pure_shells, lmax_ao);
      const auto eps = std::numeric_limits<double>::epsilon();
      Engine eng_3c(Operator::coulomb_opop, nprim, 2, 0, eps,
                    libint2::default_params(Operator::coulomb_opop),
                    BraKet::xs_xx);
      Engine eng_4c(Operator::coulomb_opop, nprim, 2, 0);
      for (size_t a = 0; a != pure.size(); ++a)
        for (size_t b = 0; b != pure.size(); ++b) {
          const std::vector<Shell> ps{fit, pure[a], pure[b]};
          const auto val =
              rkb_copy_results(eng_3c.compute(fit, pure[a], pure[b]), ps);
          const auto ref = rkb_copy_results(
              eng_4c.compute(fit, Shell::unit(), pure[a], pure[b]), ps);
          std::ostringstream oss;
          oss << "(d|" << a << b << ") 3c vs 4c";
          rkb_require_equal(val, ref, oss.str());
        }
    }
  }
#endif
#endif
#if LIBINT2_TASK_EXISTS_3coulomb_op
  SECTION("coulomb_op 3-center") {
    check_3center(Operator::coulomb_op, LIBINT2_MAX_AM_3coulomb_op,
                  LIBINT2_MAX_AM_default);
  }
#endif
}

TEST_CASE("Erfx_Coulomb integrals", "[engine][2-body]") {
#if defined(LIBINT2_SUPPORT_ERI)
  // pseudorandom s shells
  std::vector<Shell> obs{
      Shell{{1.0}, {{0, true, {1.0}}}, {{0.0, 0.0, 0.0}}},
      Shell{{2.0}, {{0, true, {1.0}}}, {{1.0, 1.0, 1.0}}},
      Shell{{3.0}, {{0, true, {1.0}}}, {{2.0, 2.0, 2.0}}},
      Shell{{4.0}, {{0, true, {1.0}}}, {{3.0, 3.0, 3.0}}},
  };
  const auto max_nprim = libint2::max_nprim(obs);
  const auto max_l = libint2::max_l(obs);

  if (LIBINT2_MAX_AM_eri >= max_l) {
    const double omega = 1.1;
    std::map<int, Operator> int_to_op = {{0, Operator::erf_coulomb},
                                         {1, Operator::erfc_coulomb},
                                         {2, Operator::erfx_coulomb}};
    for (int k = 0; k <= 2; ++k) {
      const auto op = int_to_op[k];
      auto engine = Engine(op, max_nprim, max_l);
      if (op == Operator::erf_coulomb || op == Operator::erfc_coulomb) {
        engine.set_params(omega);
      } else {  // Operator::erfx_coulomb
        engine.set_params(std::array<double, 3>{omega, 2., 3.});
      }
      const auto &results = engine.results();

      engine.compute(obs[0], obs[1], obs[2], obs[3]);
      REQUIRE(results[0] != nullptr);
      switch (k) {
        /* VALIDATION WOLFRAM CODE:
 (* Integral of Coulomb kernel damped by (\[Lambda] Erf[\[Omega] r] + \
 \[Sigma] Erfc[\[Omega] r]), over unit-normalized s functions, \
 see Eq 52 in DOI 10.1039/b605188j *)
 F0[T_] := If[T == 0, 1, Sqrt[\[Pi]/T]*Erf[Sqrt[T]]/2];
 sN[a_] := ((2 a)/\[Pi])^(3/4);
 VVeeErfx[\[Alpha]1_, A1_List, \[Alpha]2_, A2_List, \[Beta]1_,
   B1_List, \[Beta]2_, B2_List, \[Omega]_, \[Lambda]_, \[Sigma]_] :=
  Module[{\[Gamma]1, \[Gamma]2, P1, P2, K1, K2, T, result, \[Rho]},
   \[Gamma]1 = \[Alpha]1 + \[Beta]1;
   \[Gamma]2 = \[Alpha]2 + \[Beta]2;
   P1 = (\[Alpha]1*A1 + \[Beta]1*B1)/\[Gamma]1;
   P2 = (\[Alpha]2*A2 + \[Beta]2*B2)/\[Gamma]2;
   K1 = Exp[-\[Alpha]1*\[Beta]1*(Norm[A1 - B1]^2)/\[Gamma]1];
   K2 = Exp[-\[Alpha]2*\[Beta]2*(Norm[A2 - B2]^2)/\[Gamma]2];
   T = Norm[P1 - P2]^2*\[Gamma]1*\[Gamma]2/(\[Gamma]1 + \[Gamma]2);
   \[Rho] = \[Gamma]1 \[Gamma]2/(\[Gamma]1 + \[Gamma]2);
   result = (\[Pi]/(\[Gamma]1 + \[Gamma]2))^(3/
        2) K1 K2  (2 \[Pi]/\[Rho]) (\[Sigma] F0[
         T] - (\[Sigma] - \[Lambda]) \[Omega]/
        Sqrt[\[Omega]^2 + \[Rho]] F0[\[Omega]^2/(\[Omega]^2 + \[Rho])
           T]) sN[\[Alpha]1] sN[\[Alpha]2] sN[\[Beta]1] sN[\[Beta]2];
   Return[result];
   ];
 Print[CForm[
  N[VVeeErfx[1, {0, 0, 0},  3, {2, 2, 2}, 2, {1, 1, 1}, 4, {3, 3, 3},
    1.1, 1, 0], 20]]]
 Print[CForm[
  N[VVeeErfx[1, {0, 0, 0},  3, {2, 2, 2}, 2, {1, 1, 1}, 4, {3, 3, 3},
    1.1, 0, 1], 20]]]
 Print[CForm[
  N[VVeeErfx[1, {0, 0, 0},  3, {2, 2, 2}, 2, {1, 1, 1}, 4, {3, 3, 3},
    1.1, 2, 3], 20]]]
         */
        case 0:
          CHECK(results[0][0] == Approx(0.00021597118358999701).margin(1e-13));
          break;
        case 1:
          CHECK(results[0][0] == Approx(9.399863466997701e-9).margin(1e-13));
          break;
        case 2:
          CHECK(results[0][0] == Approx(0.0004319705667703951).margin(1e-13));
          break;
      }
    }
  }
#endif  // LIBINT2_SUPPORT_ERI
}
