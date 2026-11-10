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

#include <libint2/cgshell_ordering.h>
#include <libint2/shgshell_ordering.h>

#include <array>
#include <vector>

#include "catch.hpp"
#include "fixture.h"

TEST_CASE_METHOD(libint2::unit::DefaultFixture, "solid harmonics", "[engine]") {
  std::vector<Shell> obs{
      Shell{{1.0, 3.0}, {{1, true, {1.0, 0.3}}}, {{0.0, 0.0, 0.0}}},
      Shell{{1.0, 3.0},
            {{1, false, {1.0, 0.3}}},
            {{0.0, 0.0, 0.0}}},  // same as 0, but Cartesian
      Shell{{2.0, 5.0}, {{2, true, {1.0, 0.2}}}, {{1.0, 1.0, 1.0}}}};

  // test that integrals with p solid harmonics are correctly computed
  // N.B. since by default BasisSet makes p's cartesian this is rarely (if at
  // all) tested in production
  SECTION("solid harmonic p") {
#if defined(LIBINT2_SUPPORT_ONEBODY)
    const auto lmax = std::min(2, LIBINT2_MAX_AM_elecpot);
    if (lmax >= 2) {
      auto engine = Engine(Operator::nuclear, 2, lmax);
      engine.set_params(make_point_charges(atoms));

      engine.compute(obs[0], obs[2]);
      auto ints02 = engine.results()[0];
      std::vector<double> PD_ints(ints02, ints02 + 15);
      engine.compute(obs[1], obs[2]);
      auto ints12 = engine.results()[0];
      std::vector<double> pD_ints(ints12, ints12 + 15);

      int ix, iy, iz;
      int xyz = 0;
      LIBINT2_FOR_CART(ix, iy, iz, 1)
      int lm;
      if (ix == 1) {
        lm = libint2::INT_SOLIDHARMINDEX(libint2::solid_harmonics_ordering(), 1,
                                         1);
      } else if (iy == 1) {
        lm = libint2::INT_SOLIDHARMINDEX(libint2::solid_harmonics_ordering(), 1,
                                         -1);
      } else if (iz == 1) {
        lm = libint2::INT_SOLIDHARMINDEX(libint2::solid_harmonics_ordering(), 1,
                                         0);
      }

      for (int f2 = 0; f2 != 5; ++f2) {
        auto cart_f1f2 = xyz * 5 + f2;
        auto sph_f1f2 = lm * 5 + f2;
        REQUIRE(PD_ints[sph_f1f2] == Approx(pD_ints[cart_f1f2]));
      }

      ++xyz;
      LIBINT2_END_FOR_CART
    }

#endif  // LIBINT2_SUPPORT_ONEBODY
  }     // SECTION("solid harmonic p")
}

// the loop macros must work outside of namespace libint2 (issue #376), hence
// these helpers are at global scope
namespace {
std::vector<std::array<int, 3>> cart_loop(int am) {
  std::vector<std::array<int, 3>> result;
  int i, j, k;
  LIBINT2_FOR_CART(i, j, k, am)
  result.push_back({i, j, k});
  LIBINT2_END_FOR_CART
  return result;
}
std::vector<std::array<int, 3>> cart_molden_loop(int am) {
  std::vector<std::array<int, 3>> result;
  int i, j, k;
  LIBINT2_FOR_CART_MOLDEN(i, j, k, am)
  result.push_back({i, j, k});
  LIBINT2_END_FOR_CART_MOLDEN
  return result;
}
std::vector<int> solidharm_loop(int l) {
  std::vector<int> result;
  int m;
  LIBINT2_FOR_SOLIDHARM(l, m)
  result.push_back(m);
  LIBINT2_END_FOR_SOLIDHARM
  LIBINT2_FOR_SOLIDHARM_MOLDEN(l, m)
  result.push_back(m);
  LIBINT2_END_FOR_SOLIDHARM_MOLDEN
  return result;
}

// the deprecated unprefixed macros must behave exactly like their
// replacements; silence the (intended) deprecation warnings
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#elif defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#endif
std::vector<std::array<int, 3>> deprecated_cart_loop(int am) {
  std::vector<std::array<int, 3>> result;
  int i, j, k;
  FOR_CART(i, j, k, am)
  result.push_back({i, j, k});
  END_FOR_CART
  return result;
}
std::vector<std::array<int, 3>> deprecated_cart_molden_loop(int am) {
  std::vector<std::array<int, 3>> result;
  int i, j, k;
  FOR_CART_MOLDEN(i, j, k, am)
  result.push_back({i, j, k});
  END_FOR_CART_MOLDEN
  return result;
}
std::vector<int> deprecated_solidharm_loop(int l) {
  std::vector<int> result;
  int m;
  FOR_SOLIDHARM(l, m)
  result.push_back(m);
  END_FOR_SOLIDHARM
  FOR_SOLIDHARM_MOLDEN(l, m)
  result.push_back(m);
  END_FOR_SOLIDHARM_MOLDEN
  return result;
}
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#elif defined(_MSC_VER)
#pragma warning(pop)
#endif
}  // namespace

TEST_CASE("shell loop macros", "[shell-order]") {
  for (int l = 0; l <= 4; ++l) {
    const auto cart = cart_loop(l);
    REQUIRE(cart.size() == static_cast<std::size_t>(libint2::INT_NCART(l)));
    for (std::size_t f = 0; f != cart.size(); ++f) {
      REQUIRE(cart[f][0] + cart[f][1] + cart[f][2] == l);
      REQUIRE(libint2::INT_CARTINDEX(l, cart[f][0], cart[f][1]) ==
              static_cast<int>(f));
    }
    const auto cart_molden = cart_molden_loop(l);
    for (std::size_t f = 0; f != cart_molden.size(); ++f)
      REQUIRE(libint2::INT_CARTINDEX_MOLDEN(l, cart_molden[f][0],
                                            cart_molden[f][1]) ==
              static_cast<int>(f));
    REQUIRE(solidharm_loop(l).size() ==
            static_cast<std::size_t>(2 * (2 * l + 1)));

    REQUIRE(deprecated_cart_loop(l) == cart);
    REQUIRE(deprecated_cart_molden_loop(l) == cart_molden);
    REQUIRE(deprecated_solidharm_loop(l) == solidharm_loop(l));
  }
}
