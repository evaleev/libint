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

#include "catch.hpp"
#include "fixture.h"

#if defined(NO_LIBINT_COMPILER_CODE)
#include "../eri/eri.h"
#else
#include <test_eri/eri.h>
#endif

typedef unsigned int uint;

TEST_CASE("Slater/Yukawa integrals", "[engine][2-body]") {
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
}

// see https://github.com/evaleev/libint/issues/216
TEST_CASE_METHOD(libint2::unit::DefaultFixture, "bra-ket permutation",
                 "[engine][2-body]") {
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
}

TEST_CASE("eri geometric derivatives", "[engine][2-body]") {
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
}

TEST_CASE("RKB Coulomb integrals", "[engine][2-body]") {
  std::vector<Shell> obs{// pseudorandom s
                         Shell{{1.0}, {{0, false, {1.0}}}, {{0.0, 0.0, 0.0}}},
                         // pseudorandom p
                         Shell{{2.0}, {{1, false, {1.0}}}, {{1.0, 1.0, 1.0}}}};

  const auto max_nprim = libint2::max_nprim(obs);
  const auto max_l = libint2::max_l(obs);
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

  // deriv indices for (SS|SS)
  // 0th component
  der_idx xxxx = {1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0};
  der_idx yyxx = {0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0};
  der_idx zzxx = {0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0};
  der_idx yxyx = {0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0};
  der_idx xyyx = {1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0};
  der_idx yxxy = {0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0};
  der_idx xyxy = {1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0};
  der_idx xxyy = {1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0};
  der_idx yyyy = {0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0};
  der_idx zzyy = {0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0};
  der_idx xxzz = {1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1};
  der_idx yyzz = {0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1};
  der_idx zzzz = {0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1};

  // x-component
  der_idx zxzx = {0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0};
  der_idx xzzx = {1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0};
  der_idx zyzy = {0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0};
  der_idx yzzy = {0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0};
  der_idx zxxz = {0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1};
  der_idx xzxz = {1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1};
  der_idx zyyz = {0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1};
  der_idx yzyz = {0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1};

  // y-component
  der_idx zyzx = {0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0};
  der_idx yzzx = {0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0};
  der_idx zxzy = {0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0};
  der_idx xzzy = {1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0};
  der_idx zyxz = {0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1};
  der_idx yzxz = {0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1};
  der_idx zxyz = {0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1};
  der_idx xzyz = {1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1};

  // z-component
  der_idx yxxx = {0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0};
  der_idx xyxx = {1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0};
  der_idx xxyx = {1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0};
  der_idx yyyx = {0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0};
  der_idx zzyx = {0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0};
  der_idx xxxy = {1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0};
  der_idx yyxy = {0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0};
  der_idx zzxy = {0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0};
  der_idx yxyy = {0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0};
  der_idx xyyy = {1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0};
  der_idx yxzz = {0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1};
  der_idx xyzz = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1};

  SECTION("Coulombσpσp and σpσpCoulombσpσp") {
    Engine engine_llss, engine_ssss;
    try {
      engine_llss = Engine(Operator::coulomb_opop, max_nprim, max_l, 0);
      engine_ssss = Engine(Operator::opop_coulomb_opop, max_nprim, max_l, 0);
      // TODO: need another unit test for derivatives of RKB ERIs
    } catch (
        Engine::lmax_exceeded &) {  // skip the test if lmax exceeded or libint2
                                    // not configured with RKB support
      return;
    }

    const auto nshell = obs.size();
    for (int s0 = 0; s0 != nshell; ++s0) {
      for (int s1 = 0; s1 != nshell; ++s1) {
        for (int s2 = 0; s2 != nshell; ++s2) {
          for (int s3 = 0; s3 != nshell; ++s3) {
            const auto &results_llss =
                engine_llss.compute(obs[s0], obs[s1], obs[s2], obs[s3]);
            const auto &results_ssss =
                engine_ssss.compute(obs[s0], obs[s1], obs[s2], obs[s3]);
            assert(results_llss.size() ==
                   4);  // 4 buffers for single-spin quaternion components

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

            // (SS|SS) has 16 components (two independent spin spaces)
            for (auto comp = 0; comp < 16; ++comp) {
              auto abs_err_ssss =
                  abs(ref_opop_coulomb_opop[comp] - results_ssss[comp][ijkl]);
              auto rel_abs_err_ssss =
                  abs(abs_err_ssss / ref_opop_coulomb_opop[comp]);

              bool ssss_not_ok =
                  rel_abs_err_ssss > RELATIVE_DEVIATION_THRESHOLD &&
                  abs_err_ssss > ABSOLUTE_DEVIATION_THRESHOLD;

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

  SECTION("op_coulomb_op") {
    Engine engine_opCop;
    try {
      engine_opCop = Engine(Operator::op_coulomb_op, max_nprim, max_l, 0);
    } catch (Engine::lmax_exceeded &) {
      return;
    }

    const auto nshell = obs.size();
    for (int s0 = 0; s0 != nshell; ++s0) {
      for (int s1 = 0; s1 != nshell; ++s1) {
        for (int s2 = 0; s2 != nshell; ++s2) {
          for (int s3 = 0; s3 != nshell; ++s3) {
            const auto &results =
                engine_opCop.compute(obs[s0], obs[s1], obs[s2], obs[s3]);
            assert(results.size() == 9);

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
            // op_coulomb_op (must match opCoulombop_Descr::Component order).
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
}

TEST_CASE("Erfx_Coulomb integrals", "[engine][2-body]") {
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
}
