/* Independent truncated-Coulomb tests. LGPL-3.0-or-later; see COPYING.LESSER.
 */
#include <libint2/boys.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <libint2.hpp>
#include <limits>
#include <vector>

#include "catch.hpp"
#include "truncated-coulomb-reference.h"

TEST_CASE("Truncated Coulomb auxiliary integrals",
          "[core-ints][truncated-coulomb]") {
  const double eps = std::numeric_limits<double>::epsilon();
  libint2::os_core_ints::truncated_coulomb_gm_eval<double> eval(64, eps);
  std::array<double, 65> values;
  for (const auto& ref : truncated_coulomb_reference) {
    // Check that the result is independent of the requested prefix, including
    // the high-order region where direct upward recurrence loses accuracy.
    for (int mmax : {0, 1, 2, 8, 9, 21, 40, 64}) {
      eval(values.data(), 1., ref.t, mmax, ref.r);
      for (int m = 0; m <= mmax; ++m) {
        INFO("R=" << ref.r << " T=" << ref.t << " mmax=" << mmax << " m=" << m);
        CHECK(std::isfinite(values[m]));
        CHECK(std::abs(values[m] - ref.g[m]) <= eps);
      }
      // Equal dimensionless R and T at a different physical scale.
      eval(values.data(), 4., ref.t, mmax, ref.r / 2);
      for (int m = 0; m <= mmax; ++m)
        CHECK(std::abs(values[m] - ref.g[m]) <= eps);
    }
  }
  for (double t : {0., 0.1, 10., 100., 1e5}) {
    eval(values.data(), 1., t, 64, 0.);
    for (auto value : values) CHECK(value == 0.);
    eval(values.data(), 1., t, 64, std::numeric_limits<double>::infinity());
    std::array<long double, 65> ref;
    libint2::FmEval_Reference2<long double>::eval(ref.data(), t, 64);
    for (int m = 0; m <= 64; ++m) CHECK(std::abs(values[m] - ref[m]) <= eps);
  }
  CHECK_THROWS_AS(eval(values.data(), 1., 0., 1, -1.), std::invalid_argument);
  CHECK_THROWS_AS(
      eval(values.data(), 1., 0., 1, std::numeric_limits<double>::quiet_NaN()),
      std::invalid_argument);
  CHECK_THROWS_AS(eval(values.data(), 0., 0., 1, 1.), std::invalid_argument);
}

#ifdef LIBINT_INCLUDE_ERI
TEST_CASE("Truncated Coulomb Engine limits and derivatives",
          "[engine][2-body][truncated-coulomb]") {
  using namespace libint2;
  const auto op = Operator::truncated_coulomb;
  REQUIRE(rank(op) == 2);
  const double infinity = std::numeric_limits<double>::infinity();
  Shell a{{0.8}, {{0, false, {1.0}}}, {{0.1, -0.2, 0.3}}};
  Shell b{{1.1}, {{0, false, {1.0}}}, {{-0.4, 0.3, 0.2}}};
  Shell c{{0.6}, {{0, false, {1.0}}}, {{0.5, 0.1, -0.1}}};
  Shell d{{1.3}, {{0, false, {1.0}}}, {{-0.1, -0.3, 0.4}}};

  std::vector<std::pair<BraKet, int>> brakets{
      {BraKet::xx_xx, LIBINT_INCLUDE_ERI}};
#ifdef LIBINT_INCLUDE_ERI3
  brakets.emplace_back(BraKet::xs_xx, LIBINT_INCLUDE_ERI3);
#endif
#ifdef LIBINT_INCLUDE_ERI2
  brakets.emplace_back(BraKet::xs_xs, LIBINT_INCLUDE_ERI2);
#endif
  for (const auto& entry : brakets) {
    const auto braket = entry.first;
    for (int deriv = 0; deriv <= std::min(2, entry.second); ++deriv) {
      Engine trunc(op, 1, 1, deriv);
      Engine coul(Operator::coulomb, 1, 1, deriv);
      trunc.set(braket).set_params(infinity);
      coul.set(braket);
      auto compute = [&](Engine& e) {
        if (braket == BraKet::xx_xx)
          e.compute(a, b, c, d);
        else if (braket == BraKet::xs_xx)
          e.compute(a, c, d);
        else
          e.compute(a, c);
      };
      compute(trunc);
      compute(coul);
      REQUIRE(trunc.results().size() == coul.results().size());
      for (std::size_t k = 0; k < trunc.results().size(); ++k)
        CHECK(trunc.results()[k][0] ==
              Approx(coul.results()[k][0]).margin(2e-14));
      trunc.set_params(0.);
      compute(trunc);
      if (trunc.results()[0])
        for (auto target : trunc.results()) CHECK(target[0] == 0.);
    }
  }

  Engine e(op, 1, 0, 0);
  CHECK_THROWS_AS(e.set_params(-1.), std::invalid_argument);
  e.set_params(1.2);
  e.compute(a, b, c, d);
  const auto value = e.results()[0][0];
  e.compute(c, d, a, b);
  CHECK(e.results()[0][0] == Approx(value).margin(2e-15));
  e.compute(b, a, d, c);
  CHECK(e.results()[0][0] == Approx(value).margin(2e-15));

#if LIBINT_INCLUDE_ERI >= 1
  Engine grad(op, 1, 0, 1);
  grad.set_params(1.2);
  grad.compute(a, b, c, d);
  std::vector<double> analytic;
  for (auto target : grad.results()) analytic.push_back(target[0]);
  const double step = 1e-4;
  for (int center = 0; center < 4; ++center) {
    for (int xyz = 0; xyz < 3; ++xyz) {
      std::array<Shell, 4> shells{{a, b, c, d}};
      auto sample = [&](double displacement) {
        shells[center].O[xyz] += displacement;
        e.compute(shells[0], shells[1], shells[2], shells[3]);
        const auto v = e.results()[0][0];
        shells[center].O[xyz] -= displacement;
        return v;
      };
      const auto fd = (sample(-2 * step) - 8 * sample(-step) +
                       8 * sample(step) - sample(2 * step)) /
                      (12 * step);
      CHECK(analytic[3 * center + xyz] == Approx(fd).margin(3e-11));
    }
  }
#endif
#if LIBINT_INCLUDE_ERI >= 2
  Engine hessian(op, 1, 0, 2);
  hessian.set_params(1.2);
  hessian.compute(a, b, c, d);
  const double hxx = hessian.results()[0][0];
  const double h = 1e-3;
  auto displaced = [&](double dx) {
    Shell moved = a;
    moved.O[0] += dx;
    e.compute(moved, b, c, d);
    return e.results()[0][0];
  };
  const auto fd2 = (-displaced(2 * h) + 16 * displaced(h) - 30 * displaced(0) +
                    16 * displaced(-h) - displaced(-2 * h)) /
                   (12 * h * h);
  CHECK(hxx == Approx(fd2).margin(2e-9));
#endif
}

TEST_CASE("Truncated Coulomb contractions and reduced centers",
          "[engine][2-body][truncated-coulomb]") {
  using namespace libint2;
  const auto op = Operator::truncated_coulomb;
  const double inf = std::numeric_limits<double>::infinity();
  for (bool pure : {false, true}) {
    Shell a{{0.8, 1.7}, {{1, pure, {0.6, 0.4}}}, {{0.1, -0.2, 0.3}}};
    Shell b{{1.1}, {{0, pure, {1.0}}}, {{-0.4, 0.3, 0.2}}};
    Shell c{{0.6, 1.3}, {{1, pure, {0.7, 0.3}}}, {{0.5, 0.1, -0.1}}};
    Shell d{{1.3}, {{0, pure, {1.0}}}, {{-0.1, -0.3, 0.4}}};
    Engine e(op, 2, 0), coul(Operator::coulomb, 2, 1);
    // Exercise resizing the radial evaluator after constructing an s engine.
    e.set_max_l(1).set_params(inf);
    e.compute(a, b, c, d);
    coul.compute(a, b, c, d);
    const auto size = a.size() * b.size() * c.size() * d.size();
    for (std::size_t i = 0; i < size; ++i)
      CHECK(e.results()[0][i] == Approx(coul.results()[0][i]).margin(3e-14));
    e.set_params(1.2);
    e.compute(a, b, c, d);
    std::vector<double> ref(e.results()[0], e.results()[0] + size);
    e.compute(c, d, a, b);
    for (std::size_t i = 0; i < a.size(); ++i)
      for (std::size_t j = 0; j < c.size(); ++j)
        CHECK(e.results()[0][j * a.size() + i] ==
              Approx(ref[i * c.size() + j]).margin(2e-14));

#ifdef LIBINT_INCLUDE_ERI3
    Engine three(op, 2, 1);
    three.set(BraKet::xs_xx).set_params(1.2);
    e.compute(a, Shell::unit(), c, d);
    ref.assign(e.results()[0], e.results()[0] + a.size() * c.size() * d.size());
    three.compute(a, c, d);
    for (std::size_t i = 0; i < ref.size(); ++i)
      CHECK(three.results()[0][i] == Approx(ref[i]).margin(2e-14));
#endif
#ifdef LIBINT_INCLUDE_ERI2
    Engine two(op, 2, 1);
    two.set(BraKet::xs_xs).set_params(1.2);
    e.compute(a, Shell::unit(), c, Shell::unit());
    ref.assign(e.results()[0], e.results()[0] + a.size() * c.size());
    two.compute(a, c);
    for (std::size_t i = 0; i < ref.size(); ++i)
      CHECK(two.results()[0][i] == Approx(ref[i]).margin(2e-14));
    // Increase center count after creating a two-center core evaluator.
    two.set(BraKet::xx_xx);
    two.compute(a, b, c, d);
    e.compute(a, b, c, d);
    for (std::size_t i = 0; i < size; ++i)
      CHECK(two.results()[0][i] == Approx(e.results()[0][i]).margin(2e-14));
#endif
  }
}

#endif  // LIBINT_INCLUDE_ERI
