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

#include <libint2/solidharmonics.h>

#include "catch.hpp"
#include "fixture.h"

TEST_CASE("SolidHarmonicsCoefficients ctor", "[shell]") {
  using libint2::solidharmonics::SolidHarmonicsCoefficients;
  REQUIRE_NOTHROW(SolidHarmonicsCoefficients<long double>::instance(12));
  auto sh12 = SolidHarmonicsCoefficients<long double>::instance(12);
  CHECK_NOTHROW(sh12.coeff(12, 0, 0, 0, 12));
  CHECK(sh12.coeff(12, 0, 0, 0, 12) == Approx(1));
  CHECK_NOTHROW(
      SolidHarmonicsCoefficients<long double>::make_instance(25).coeff(25, 0, 0,
                                                                       0, 25));
  CHECK(SolidHarmonicsCoefficients<long double>::make_instance(25).coeff(
            25, 0, 0, 0, 25) == Approx(1));
}

TEST_CASE("Shell ctor", "[shell]") {
  REQUIRE_NOTHROW(Shell{});
  auto s0 = Shell{};
  REQUIRE(s0.alpha.empty());
  REQUIRE(s0.contr.empty());
  REQUIRE(s0.max_ln_coeff.empty());

  REQUIRE_NOTHROW(Shell{{1}, {{2, false, {1}}}, {{0, 0, 0}}});
  auto s1 = Shell{{1}, {{2, false, {1}}}, {{0, 0, 0}}};
  REQUIRE(s1.alpha == libint2::svector<double>{1});
  REQUIRE(s1.contr.size() == 1);
  REQUIRE(s1.contr[0].l == 2);
  REQUIRE(!s1.contr[0].pure);
  REQUIRE(s1.contr[0].coeff.size() == 1);
  REQUIRE(s1.contr[0].coeff[0] ==
          Approx(1.64592278064949));  // (2./\[Pi])^(3/4) 2^l
                                      // \[Alpha]^((2l+3)/4) / Sqrt[(2l-1)!!]
  REQUIRE(s1.O == std::array<double, 3>{0, 0, 0});
}

TEST_CASE("Engine ctor", "[engine]") {
  REQUIRE_NOTHROW(Engine{});
  auto a = Engine{};
}

TEST_CASE("Engine::set", "[engine]") {
  REQUIRE_THROWS_AS(Engine{}.set(Operator::overlap),
                    Engine::using_default_initialized);
  REQUIRE_THROWS_AS(Engine{}.set(BraKet::x_x),
                    Engine::using_default_initialized);
  REQUIRE_NOTHROW(Engine{}.set(CartesianShellNormalization::uniform));
  REQUIRE_NOTHROW(Engine{}.prescale_by(1.5));
  REQUIRE_NOTHROW(Engine(Operator::overlap, 1, 0)
                      .set(CartesianShellNormalization::uniform)
                      .set_precision(1e-20)
                      .set(Operator::overlap)
                      .prescale_by(1.3)
                      .set(BraKet::x_x));
}

template <typename RealBuf>
void check_uniform(int l, RealBuf&& S) {
  const auto n = (l + 1) * (l + 2) / 2;
  for (int i = 0; i != n; ++i) {
    REQUIRE(S[i * n + i] == Approx(1.));
  }
};

template <typename RealBuf>
void check_std_2(RealBuf&& S) {
  REQUIRE(S[0] == Approx(1.));
  REQUIRE(S[7] == Approx(1. / 3));
  REQUIRE(S[14] == Approx(1. / 3));
  REQUIRE(S[21] == Approx(1.));
  REQUIRE(S[28] == Approx(1. / 3));
  REQUIRE(S[35] == Approx(1.));
};

template <typename RealBuf>
void check_std_3(RealBuf&& S) {
  REQUIRE(S[0] == Approx(1.));
  REQUIRE(S[11] == Approx(1. / 5));
  REQUIRE(S[22] == Approx(1. / 5));
  REQUIRE(S[33] == Approx(1. / 5));
  REQUIRE(S[44] == Approx(1. / 15));
  REQUIRE(S[55] == Approx(1. / 5));
  REQUIRE(S[66] == Approx(1.));
  REQUIRE(S[77] == Approx(1. / 5));
  REQUIRE(S[88] == Approx(1. / 5));
  REQUIRE(S[99] == Approx(1.));
};

TEST_CASE("cartesian uniform normalization", "[engine][conventions]") {
#if defined(LIBINT2_SUPPORT_ONEBODY) && defined(LIBINT2_SUPPORT_ERI)
  if (LIBINT_CGSHELL_ORDERING != LIBINT_CGSHELL_ORDERING_STANDARD) return;

  std::vector<Shell> obs{Shell{{1.0}, {{2, false, {1.0}}}, {{0.0, 0.0, 0.0}}},
                         Shell{{1.0}, {{3, false, {1.0}}}, {{0.0, 0.0, 0.0}}}};
  {
    const auto lmax = std::min(3, LIBINT2_MAX_AM_overlap);
    if (lmax >= 2) {
      auto engine = Engine(Operator::overlap, 1, lmax);
      engine.compute(obs[0], obs[0]);
      check_std_2(engine.results()[0]);
      engine.set(CartesianShellNormalization::uniform).compute(obs[0], obs[0]);
      check_uniform(2, engine.results()[0]);

      if (lmax >= 3) {
        engine.set(CartesianShellNormalization::standard)
            .compute(obs[1], obs[1]);
        check_std_3(engine.results()[0]);
        engine.set(CartesianShellNormalization::uniform)
            .compute(obs[1], obs[1]);
        check_uniform(3, engine.results()[0]);
      }
    }
  }
  {
    const auto lmax = std::min(3, LIBINT2_MAX_AM_eri);
    if (lmax >= 2) {
      auto engine = Engine(Operator::delta, 1, lmax);
      engine.compute(Shell::unit(), obs[0], obs[0], Shell::unit());
      check_std_2(engine.results()[0]);
      engine.set(CartesianShellNormalization::uniform)
          .compute(obs[0], Shell::unit(), Shell::unit(), obs[0]);
      check_uniform(2, engine.results()[0]);

      if (lmax >= 3) {
        engine.set(CartesianShellNormalization::standard)
            .compute(Shell::unit(), obs[1], Shell::unit(), obs[1]);
        check_std_3(engine.results()[0]);
        engine.set(CartesianShellNormalization::uniform)
            .compute(obs[1], Shell::unit(), obs[1], Shell::unit());
        check_uniform(3, engine.results()[0]);
      }
    }
  }

#endif  // LIBINT2_SUPPORT_ONEBODY
}
