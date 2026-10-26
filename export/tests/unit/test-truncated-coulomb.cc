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
#include "truncated-coulomb-vrr-reference.h"

namespace {
template <typename Evaluator>
void check_truncated_coulomb_radial(Evaluator& eval) {
  const double eps = std::numeric_limits<double>::epsilon();
  std::array<double, 65> values;
  for (const auto& ref : truncated_coulomb_reference) {
    // Check that the result is independent of the requested prefix, including
    // the high-order region where direct upward recurrence loses accuracy.
    for (int mmax : {0, 1, 2, 8, 9, 21, 40, 64}) {
      eval.eval(values.data(), 1., ref.t, mmax, ref.r);
      for (int m = 0; m <= mmax; ++m) {
        INFO("R=" << ref.r << " T=" << ref.t << " mmax=" << mmax << " m=" << m);
        CHECK(std::isfinite(values[m]));
        CHECK(std::abs(values[m] - ref.g[m]) <= eps);
      }
      // Equal dimensionless R and T at a different physical scale.
      eval.eval(values.data(), 4., ref.t, mmax, ref.r / 2);
      for (int m = 0; m <= mmax; ++m)
        CHECK(std::abs(values[m] - ref.g[m]) <= eps);
    }
  }
  for (double t : {0., 0.1, 10., 100., 1e5}) {
    eval.eval(values.data(), 1., t, 64, 0.);
    for (auto value : values) CHECK(value == 0.);
    eval.eval(values.data(), 1., t, 64,
              std::numeric_limits<double>::infinity());
    std::array<long double, 65> ref;
    libint2::FmEval_Reference2<long double>::eval(ref.data(), t, 64);
    for (int m = 0; m <= 64; ++m) CHECK(std::abs(values[m] - ref[m]) <= eps);
  }
  CHECK_THROWS_AS(eval.eval(values.data(), 1., 0., 1, -1.),
                  std::invalid_argument);
  CHECK_THROWS_AS(eval.eval(values.data(), 1., 0., 1,
                            std::numeric_limits<double>::quiet_NaN()),
                  std::invalid_argument);
  CHECK_THROWS_AS(eval.eval(values.data(), 0., 0., 1, 1.),
                  std::invalid_argument);
}
}  // namespace

TEST_CASE("Truncated Coulomb auxiliary integrals",
          "[core-ints][truncated-coulomb]") {
  const double eps = std::numeric_limits<double>::epsilon();
  libint2::os_core_ints::truncated_coulomb_gm_eval<double> tc(192, eps);
  check_truncated_coulomb_radial(tc);

  // Independent defining-erf references, generated at 100 decimal digits and
  // checked at 140 digits. These inputs exposed cancellation at m=40 and a
  // parasitic high-order solution at m=192; the criterion is absolute error.
  std::array<double, 193> values;
  decltype(tc)::scratch_type scratch(192);
  tc.eval(values.data(), 1., 15., 40, 1., scratch);
  CHECK(std::abs(values[40] - 9.1989225488094219e-09) <= eps);
  tc.eval(values.data(), 1., 59., 192, .1, scratch);
  CHECK(std::abs(values[192] - 1.4877535827085252e-28) <= eps);
  // Reusing a larger workspace must not retain the preceding hierarchy.
  tc.eval(values.data(), 1., 15., 40, 1., scratch);
  CHECK(std::abs(values[40] - 9.1989225488094219e-09) <= eps);
  // A reserve beyond the native Boys table must still use its supported
  // orders for later, lower-order requests, even if the singleton has grown.
  std::array<double, 9> coulomb;
  libint2::FmEval_Chebyshev7<double>::instance(8)->eval(coulomb.data(), .24, 8);
  tc.eval(values.data(), 1., .24, 8, std::numeric_limits<double>::infinity(),
          scratch);
  for (int m = 0; m <= 8; ++m) CHECK(values[m] == coulomb[m]);
  tc.eval(values.data(), 1., .24, 8, 20., scratch);
  for (int m = 0; m <= 8; ++m) CHECK(std::abs(values[m] - coulomb[m]) <= eps);
}

// 0 < T < min(n + 1/2, 12) uses the double-word origin expansion: its top
// of range at a high order, and R well beyond sqrt(T), where U_k and V_k
// grow with the dominant solution and the error estimate needs its nested
// form. Defining-erf references at 50 digits; the criterion is absolute.
// Requests with T >= n + 1/2 take the other paths and must agree.
TEST_CASE("Truncated Coulomb auxiliary integrals by the origin expansion",
          "[core-ints][truncated-coulomb]") {
  const double eps = std::numeric_limits<double>::epsilon();
  struct Reference {
    double r, t;
    int n;
    std::vector<double> g;
  };
  const Reference refs[] = {
      {10.17,
       4.39,
       8,
       {0.42168486463540115, 0.046615505172546444, 0.014515465409703409,
        0.0068538266318990484, 0.0040519427338097875, 0.0027410883125107839,
        0.0020217815737101615, 0.0015811425099987123, 0.0012888847869632435}},
      {2.5,
       0.75,
       16,
       {0.78684342664308601, 0.23274225857973911, 0.10702959412734653,
        0.079677306351208520, 0.064820966424317426, 0.050486149778876075,
        0.039215987596959884, 0.032395090517522593, 0.028933570819893606,
        0.026978128157491320, 0.025178016413891936, 0.023023311904602649,
        0.020622817168460760, 0.018330621074997705, 0.016456474963152301,
        0.015131173955875089, 0.014302465367542966}},
      {4.0,
       6.0,
       24,
       {0.35648410865939795,    0.033264962713753418,   0.0057267648960145567,
        0.0031756712700991371,  0.0015358745407417276,  0.00079818953623173328,
        0.00057881633889808070, 0.00045803099002024138, 0.00034837190263657074,
        0.00027074722046168703, 0.00022641815596211149, 0.00019844919959376965,
        0.00017484643544159914, 0.00015342899987390319, 0.00013583485881713887,
        0.00012261888641286325, 0.00011258053758827510, 0.00010414278242677521,
        9.6430909906661019e-5,  8.9356055604776100e-5,  8.3142386188761031e-5,
        7.7896557283820035e-5,  7.3481447750625323e-5,  6.9627983537232975e-5,
        6.6101460978547776e-5}},
      {1.3,
       11.9,
       32,
       {0.00030356804207075164, 0.00021926625544449961, 0.00015621835795788411,
        0.00010956135496208932, 7.5444668313430262e-5,  5.0837636044396094e-5,
        3.3370582399328925e-5,  2.1204065307623918e-5,  1.2921709858899920e-5,
        7.4427271828251431e-6,  3.9508125676048421e-6,  1.8366272052554936e-6,
        6.5150563597032163e-7,  7.0405436508649165e-8,  -1.3756457694838029e-7,
        -1.3242897395803016e-7, -2.1602980738434684e-8, 1.2591410098219692e-7,
        2.6861185930436844e-7,  3.8415069602994759e-7,  4.6311306128179859e-7,
        5.0438309904018694e-7,  5.1179632360591788e-7,  4.9176966149080371e-7,
        4.5167919355463708e-7,  3.9880001985529897e-7,  3.3966139346921833e-7,
        2.7970195820779460e-7,  2.2313571606541575e-7,  1.7296021462410085e-7,
        1.3105520979003631e-7,  9.8333427916515778e-8,  7.4915618858296132e-8}},
  };
  libint2::os_core_ints::truncated_coulomb_gm_eval<double> tc(32, eps);
  decltype(tc)::scratch_type scratch(32);
  std::array<double, 33> g;
  for (const auto& ref : refs) {
    for (int n : {ref.n, static_cast<int>(ref.t - .5)}) {
      tc.eval(g.data(), 1., ref.t, n, ref.r, scratch);
      for (int m = 0; m <= n; ++m) {
        INFO("R=" << ref.r << " T=" << ref.t << " n=" << n << " m=" << m);
        CHECK(std::abs(g[m] - ref.g[m]) <= eps);
      }
    }
  }
}

// Just beyond the distance |R - sqrt(T)| at which an absolute criterion
// would switch to F_m or zero; the high orders are far below eps here, and
// G_16 even exceeds F_16 by 1e8 at T = 1e5. Defining-erf references at 50
// digits; the criterion is relative to the larger of |G_m| and F_m.
TEST_CASE("Truncated Coulomb auxiliary integrals near the Boys and zero limits",
          "[core-ints][truncated-coulomb]") {
  constexpr int mmax = 16;
  struct Reference {
    double r, t, g[mmax + 1];
  };
  const Reference refs[] = {
      {108.6,
       1e4,
       {0.0088622692545275801, 4.4311346272637901e-7, 6.6467019408956851e-11,
        1.6616754852239213e-14, 5.8158641982837245e-18, 2.617138889227676e-21,
        1.4394263890752218e-24, 9.3562715289889424e-28, 7.0172036467416507e-31,
        5.9646230997349737e-34, 5.66639194437939e-37, 5.9497115711008598e-40,
        6.8421659681667526e-43, 8.5528911216056172e-46, 1.1532116269244219e-48,
        1.7822050813808079e-51, -5.6290925324490396e-54}},
      {91.4,
       1e4,
       {2.1884352861136341e-36, 1.8957055345478377e-37, 1.6322092908399323e-38,
        1.3966766794916191e-39, 1.1876156724097307e-40, 1.0033637183153641e-41,
        8.4213758114545613e-43, 7.0208024419088905e-44, 5.8130393957925464e-45,
        4.7792858020065692e-46, 3.9011509332364499e-47, 3.1609337260052282e-48,
        2.5418476958424058e-49, 2.0281923028827299e-50, 1.6054734130845148e-51,
        1.2604767988555799e-52, 9.8129961037524071e-54}},
      {324.8277660168379,
       1e5,
       {0.0028024956081989643, 1.4012478040994822e-8, 2.1018717061492233e-13,
        5.2546792653730582e-18, 1.8391377428805704e-22, 8.2761198429625666e-27,
        4.551865913629409e-31, 2.9587128438598042e-35, 2.2190346327157234e-39,
        1.8861794840531705e-43, 1.7918586678269412e-47, 1.8844588437658945e-51,
        1.4099508446526933e-55, 1.9074151329311718e-57, -4.648509778688959e-59,
        1.138728452152227e-60, -2.7627424580754614e-62}},
      {307.62776601683794,
       1e5,
       {6.9204400160012852e-37, 1.8949573389055659e-38, 5.1553492922574712e-40,
        1.3933279360641161e-41, 3.7404711744624385e-43, 9.9727609630811688e-45,
        2.6403173634543515e-46, 6.9403414288238827e-48, 1.8110001939567021e-49,
        4.6902331769513531e-51, 1.2053967397591896e-52, 3.0735642877976119e-54,
        7.7739850384888677e-56, 1.9500309936894007e-57, 4.8499523714967992e-59,
        1.1957136764392447e-60, 2.9214574907968413e-62}},
  };
  libint2::os_core_ints::truncated_coulomb_gm_eval<double> tc(
      mmax, std::numeric_limits<double>::epsilon());
  decltype(tc)::scratch_type scratch(mmax);
  for (const auto& ref : refs) {
    double g[mmax + 1];
    long double f[mmax + 1];
    tc.eval(g, 1., ref.t, mmax, ref.r, scratch);
    libint2::FmEval_Reference2<long double>::eval(f, ref.t, mmax);
    for (int m = 0; m <= mmax; ++m) {
      INFO("R=" << ref.r << " T=" << ref.t << " m=" << m);
      const double scale =
          std::max(std::abs(ref.g[m]), static_cast<double>(f[m]));
      CHECK(std::abs(g[m] - ref.g[m]) <= 1e-12 * scale);
    }
  }
}

// (L0|00)^(0), L = 0..mmax, by the Obara-Saika VRR with P = A, from G_m.
template <std::size_t N>
std::array<double, N> vrr_l000(const double* g, double zeta, double rho,
                               double wp) {
  std::array<std::array<double, N>, N> a{};
  for (std::size_t m = 0; m < N; ++m) a[0][m] = g[m];
  for (std::size_t l = 0; l + 1 < N; ++l)
    for (std::size_t m = 0; m + l + 1 < N; ++m)
      a[l + 1][m] =
          wp * a[l][m + 1] +
          (l ? l / (2 * zeta) * (a[l - 1][m] - rho / zeta * a[l - 1][m + 1])
             : 0.);
  std::array<double, N> result;
  for (std::size_t l = 0; l < N; ++l) result[l] = a[l][0];
  return result;
}

// What must hold is the accuracy of the integrals, which neither a bound
// relative to each F_m nor an absolute one is the right measure of: far
// outside the cutoff band the high-order G_m enter integrals multiplied by
// small powers of the geometry, while near it an error of eps F_0 in G_4 can
// already exceed eps (L0|00) by 1e8 for tight primitives (the third
// reference geometry). Compare (L0|00) from evaluated and reference G_m,
// relative to the larger of the truncated and the Coulomb integral, against
// the double-precision VRR floor (~2^(L/2) eps; the tolerance leaves a factor
// of 8).
TEST_CASE("Truncated Coulomb integrals across the cutoff transition",
          "[core-ints][truncated-coulomb]") {
  constexpr int mmax = 16;
  constexpr std::size_t N = mmax + 1;
  const double eps = std::numeric_limits<double>::epsilon();
  libint2::os_core_ints::truncated_coulomb_gm_eval<double> tc(mmax, eps);
  decltype(tc)::scratch_type scratch(mmax);
  const auto boys = libint2::FmEval_Chebyshev7<double>::instance(mmax);
  for (const auto& ref : truncated_coulomb_vrr_reference) {
    double g[N], f[N];
    tc.eval(g, 1., ref.t, mmax, ref.r, scratch);
    boys->eval(f, ref.t, mmax);
    const auto computed = vrr_l000<N>(g, ref.zeta, ref.rho, ref.wp);
    const auto exact = vrr_l000<N>(ref.g, ref.zeta, ref.rho, ref.wp);
    const auto coulomb = vrr_l000<N>(f, ref.zeta, ref.rho, ref.wp);
    for (int l = 0; l <= mmax; ++l) {
      INFO("zeta=" << ref.zeta << " rho=" << ref.rho << " R=" << ref.r
                   << " T=" << ref.t << " L=" << l);
      const double scale = std::max(std::abs(exact[l]), std::abs(coulomb[l]));
      CHECK(std::abs(computed[l] - exact[l]) <=
            8 * eps * std::ldexp(1., l / 2) * scale);
    }
  }
}

#ifdef LIBINT_INCLUDE_ERI
TEST_CASE("Truncated Coulomb Engine limits and derivatives",
          "[engine][2-body][truncated-coulomb]") {
  using namespace libint2;
  for (const auto op : {Operator::truncated_coulomb}) {
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

    // Changing the operator must restore its default (zero) cutoff, rather
    // than retaining the old core parameters.
    e.set(Operator::coulomb);
    e.set(op);
    e.compute(a, b, c, d);
    CHECK((!e.results()[0] || e.results()[0][0] == 0.));
    e.set_params(1.2);
    e.compute(a, b, c, d);
    CHECK(e.results()[0][0] == value);

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
    const auto fd2 =
        (-displaced(2 * h) + 16 * displaced(h) - 30 * displaced(0) +
         16 * displaced(-h) - displaced(-2 * h)) /
        (12 * h * h);
    CHECK(hxx == Approx(fd2).margin(2e-9));
#endif
  }
}

// (aa|bb) of normalized s primitives at distance D equals
// [erf(sqrt(mu)(Rc-D)) - erf(sqrt(mu)(Rc+D)) + 2 erf(sqrt(mu) D)] / (2D),
// mu = 4ab/(2a+2b): the radial average of theta(Rc-r)/r over the Gaussian
// distribution of r12. This checks the cutoff units and its sqrt(rho)
// scaling independently of the radial evaluator. mpmath values, 50 digits.
TEST_CASE("Truncated Coulomb Engine against closed-form (ss|ss)",
          "[engine][2-body][truncated-coulomb]") {
  using namespace libint2;
  Shell a{{0.8}, {{0, false, {1.0}}}, {{0, 0, 0}}};
  Shell b{{1.3}, {{0, false, {1.0}}}, {{0, 0, 1.5}}};
  const std::pair<double, double> refs[] = {{0.5, 0.031552595422039718},
                                            {1.2, 0.2011613560258891},
                                            {1.5, 0.31017100129777589},
                                            {2.5, 0.59039957931514569},
                                            {4.0, 0.64335169847410041}};
  Engine e(Operator::truncated_coulomb, 1, 0);
  for (const auto& ref : refs) {
    INFO("Rc=" << ref.first);
    e.set_params(ref.first);
    e.compute(a, a, b, b);
    CHECK(e.results()[0][0] == Approx(ref.second).epsilon(1e-14));
  }
}

TEST_CASE("Truncated Coulomb contractions and reduced centers",
          "[engine][2-body][truncated-coulomb]") {
  using namespace libint2;
  for (const auto op : {Operator::truncated_coulomb}) {
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
      ref.assign(e.results()[0],
                 e.results()[0] + a.size() * c.size() * d.size());
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
}

#endif  // LIBINT_INCLUDE_ERI
