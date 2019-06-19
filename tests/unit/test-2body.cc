#include "catch.hpp"
#include "fixture.h"

TEST_CASE("Slater/Yukawa integrals", "[engine][2-body]") {

  std::vector<Shell> obs{
      Shell{{1.0, 3.0}, {{2, true, {1.0, 0.3}}}, {{0.0, 0.0, 0.0}}},
      Shell{{2.0, 5.0}, {{2, true, {1.0, 0.2}}}, {{1.0, 1.0, 1.0}}}};

  // 6-term GTG fit of exp(-r12) from DOI 10.1063/1.1999632
  std::vector<std::pair<double, double>> cgtg_params{{0.2209,0.3144},{1.004,0.3037},{3.622,0.1681},{12.16,0.09811},{45.87,0.06024},{254.4,0.03726}};

  for(int k=-1; k<=0; ++k) {
    using namespace libint2::unit;

    const auto lmax = std::min(3, LIBINT2_MAX_AM_eri);
    auto engine = Engine(k == 0 ? Operator::stg : Operator::yukawa, 2, lmax);
    engine.set_params(1.0);
    auto engine_ref = Engine(k == 0 ? Operator::cgtg : Operator::cgtg_x_coulomb, 2, lmax);
    engine_ref.set_params(cgtg_params);

    engine.compute(obs[0], obs[0], obs[0], obs[0]);
    engine_ref.compute(obs[0], obs[0], obs[0], obs[0]);
    {
      for (int i = 0; i != 625; ++i) {
        REQUIRE(engine.results()[0][i] == Approx(engine_ref.results()[0][i]).margin(1e-3));
      }
    }

    engine.compute(obs[0], obs[1], obs[0], obs[1]);
    engine_ref.compute(obs[0], obs[1], obs[0], obs[1]);
    {
      for (int i = 0; i != 625; ++i) {
        REQUIRE(engine.results()[0][i] == Approx(engine_ref.results()[0][i]).margin(1e-3));
      }
    }
  }

}
