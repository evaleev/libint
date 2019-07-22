#include "catch.hpp"
#include "fixture.h"

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
  REQUIRE(s1.contr[0].pure == false);
  REQUIRE(s1.contr[0].coeff.size() == 1);
  REQUIRE(s1.contr[0].coeff[0] == Approx(1.64592278064949)); // (2./\[Pi])^(3/4) 2^l \[Alpha]^((2l+3)/4) / Sqrt[(2l-1)!!]
  REQUIRE(s1.O == std::array<double,3>{0, 0, 0});
}

TEST_CASE("Engine ctor", "[engine]") {
  REQUIRE_NOTHROW(Engine{});
  auto a = Engine{};
}

TEST_CASE("Engine::set", "[engine]") {
  REQUIRE_THROWS_AS(Engine{}.set(Operator::overlap), Engine::using_default_initialized);
  REQUIRE_THROWS_AS(Engine{}.set(BraKet::x_x), Engine::using_default_initialized);
  REQUIRE_NOTHROW(Engine{}.set(CartesianShellNormalization::uniform));
  REQUIRE_NOTHROW(Engine(Operator::overlap, 1, 0).set(CartesianShellNormalization::uniform).set_precision(1e-20).set(Operator::overlap).set(BraKet::x_x));
}

template <typename RealBuf> void check_uniform(int l, RealBuf && S) {
  const auto n = (l+1)*(l+2)/2;
  for(int i=0; i!=n; ++i) {
    REQUIRE(S[i*n+i] == Approx(1.));
  }
};

template <typename RealBuf> void check_std_2(RealBuf && S) {
  REQUIRE(S[0] == Approx(1.));
  REQUIRE(S[7] == Approx(1./3));
  REQUIRE(S[14] == Approx(1./3));
  REQUIRE(S[21] == Approx(1.));
  REQUIRE(S[28] == Approx(1./3));
  REQUIRE(S[35] == Approx(1.));
};

template <typename RealBuf> void check_std_3(RealBuf && S) {
  REQUIRE(S[0] == Approx(1.));
  REQUIRE(S[11] == Approx(1./5));
  REQUIRE(S[22] == Approx(1./5));
  REQUIRE(S[33] == Approx(1./5));
  REQUIRE(S[44] == Approx(1./15));
  REQUIRE(S[55] == Approx(1./5));
  REQUIRE(S[66] == Approx(1.));
  REQUIRE(S[77] == Approx(1./5));
  REQUIRE(S[88] == Approx(1./5));
  REQUIRE(S[99] == Approx(1.));
};

TEST_CASE("cartesian uniform normalization", "[engine][conventions]") {
#if defined(LIBINT2_SUPPORT_ONEBODY) && defined(LIBINT2_SUPPORT_ERI)
  if (LIBINT_CGSHELL_ORDERING != LIBINT_CGSHELL_ORDERING_STANDARD)
    return;

  std::vector<Shell> obs{Shell{{1.0}, {{2, false, {1.0}}}, {{0.0, 0.0, 0.0}}},
                         Shell{{1.0}, {{3, false, {1.0}}}, {{0.0, 0.0, 0.0}}}};
  {
    const auto lmax = std::min(3,LIBINT2_MAX_AM_overlap);
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
    const auto lmax = std::min(3,LIBINT2_MAX_AM_eri);
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
