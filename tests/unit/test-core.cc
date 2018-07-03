#include "catch.hpp"
#include "fixture.h"

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

TEST_CASE("cartesian uniform normalization", "[engine]") {
#ifdef LIBINT2_SUPPORT_ONEBODY
  if (LIBINT_CGSHELL_ORDERING != LIBINT_CGSHELL_ORDERING_STANDARD)
    return;

  const auto lmax = std::min(3,LIBINT2_MAX_AM_overlap);
  auto engine = Engine(Operator::overlap, 1, lmax);
  std::vector<Shell> obs{Shell{{1.0}, {{2, false, {1.0}}}, {{0.0, 0.0, 0.0}}},
                         Shell{{1.0}, {{3, false, {1.0}}}, {{0.0, 0.0, 0.0}}}};
  engine.compute(obs[0], obs[0]);
  REQUIRE(engine.results()[0][0] == Approx(1.));
  REQUIRE(engine.results()[0][7] == Approx(1./3));
  REQUIRE(engine.results()[0][14] == Approx(1./3));
  REQUIRE(engine.results()[0][21] == Approx(1.));
  REQUIRE(engine.results()[0][28] == Approx(1./3));
  REQUIRE(engine.results()[0][35] == Approx(1.));
  engine.set(CartesianShellNormalization::uniform).compute(obs[0], obs[0]);
  REQUIRE(engine.results()[0][0] == Approx(1.));
  REQUIRE(engine.results()[0][7] == Approx(1.));
  REQUIRE(engine.results()[0][14] == Approx(1.));
  REQUIRE(engine.results()[0][21] == Approx(1.));
  REQUIRE(engine.results()[0][28] == Approx(1.));
  REQUIRE(engine.results()[0][35] == Approx(1.));

  if (lmax >= 3) {
    engine.set(CartesianShellNormalization::standard).compute(obs[1], obs[1]);
    REQUIRE(engine.results()[0][0] == Approx(1.));
    REQUIRE(engine.results()[0][11] == Approx(1./5));
    REQUIRE(engine.results()[0][22] == Approx(1./5));
    REQUIRE(engine.results()[0][33] == Approx(1./5));
    REQUIRE(engine.results()[0][44] == Approx(1./15));
    REQUIRE(engine.results()[0][55] == Approx(1./5));
    REQUIRE(engine.results()[0][66] == Approx(1.));
    REQUIRE(engine.results()[0][77] == Approx(1./5));
    REQUIRE(engine.results()[0][88] == Approx(1./5));
    REQUIRE(engine.results()[0][99] == Approx(1.));
    engine.set(CartesianShellNormalization::uniform).compute(obs[1], obs[1]);
    REQUIRE(engine.results()[0][0] == Approx(1.));
    REQUIRE(engine.results()[0][11] == Approx(1.));
    REQUIRE(engine.results()[0][22] == Approx(1.));
    REQUIRE(engine.results()[0][33] == Approx(1.));
    REQUIRE(engine.results()[0][44] == Approx(1.));
    REQUIRE(engine.results()[0][55] == Approx(1.));
    REQUIRE(engine.results()[0][66] == Approx(1.));
    REQUIRE(engine.results()[0][77] == Approx(1.));
    REQUIRE(engine.results()[0][88] == Approx(1.));
    REQUIRE(engine.results()[0][99] == Approx(1.));
  }
#endif  // LIBINT2_SUPPORT_ONEBODY
}
