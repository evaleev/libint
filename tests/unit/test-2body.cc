#include "catch.hpp"
#include "fixture.h"

TEST_CASE("Slater/Yukawa integrals", "[engine][2-body]") {

  std::vector<Shell> obs{
      // pseudorandom d
      Shell{{1.0, 3.0}, {{2, true, {1.0, 0.3}}}, {{0.0, 0.0, 0.0}}},
      // pseudorandom d
      Shell{{2.0, 5.0}, {{2, true, {1.0, 0.2}}}, {{1.0, 1.0, 1.0}}},
      // O 1s in cc-pVDZ
      Shell{{11720.0000000, 1759.0000000, 400.8000000, 113.7000000, 37.0300000,
             13.2700000, 5.0250000, 1.0130000},
            {{0,
              false,
              {0.0007100, 0.0054700, 0.0278370, 0.1048000, 0.2830620, 0.4487190,
               0.2709520, 0.0154580}}},
            {{2.0, 2.0, -1.0}}},
      // O 1s in STO-3G
      Shell{{130.709320000, 23.808861000, 6.443608300},
            {{0, false, {0.15432897, 0.53532814, 0.44463454}}},
            {{-1.0, -1.0, 0.0}}}};
  const auto max_nprim = libint2::max_nprim(obs);
  const auto max_l = libint2::max_l(obs);

  // 6-term GTG fit of exp(-r12) from DOI 10.1063/1.1999632
  //std::vector<std::pair<double, double>> cgtg_params{{0.2209,0.3144},{1.004,0.3037},{3.622,0.1681},{12.16,0.09811},{45.87,0.06024},{254.4,0.03726}};
  // precise fit of exp(-r12) on [0,10)
  std::vector<std::pair<double,double>> cgtg_params =
    {{0.10535330565471572,0.08616353459042002},{0.22084823136587992,0.08653979627551414},{0.3543431104992702,0.08803697599356214},{0.48305514665749105,0.09192519612306953},{0.6550035700167584,0.10079776426873248},{1.1960917050801643,0.11666110644901695},{2.269278814810891,0.14081371547404428},{5.953990617813977,0.13410255216448014},{18.31911063199608,0.0772095196191394},{66.98443868169818,0.049343985939540556},{367.24137290439205,0.03090625839896873},{5.655142311118115,-0.017659052507938647}};

  if (LIBINT2_MAX_AM_eri >= max_l) {
    for (int k = -1; k <= 0; ++k) {
      auto engine = Engine(k == 0 ? Operator::stg : Operator::yukawa, max_nprim, max_l);
      engine.set_params(1.0);
      const auto &results = engine.results();
      auto engine_ref =
          Engine(k == 0 ? Operator::cgtg : Operator::cgtg_x_coulomb, max_nprim, max_l);
      engine_ref.set_params(cgtg_params);
      const auto &results_ref = engine.results();

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
                  REQUIRE(results[0][i] ==
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
