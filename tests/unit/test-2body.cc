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
    const auto deriv_order = INCLUDE_ERI;
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
