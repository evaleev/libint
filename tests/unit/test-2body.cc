#include "catch.hpp"
#include "fixture.h"

TEST_CASE("Slater/Yukawa integrals", "[engine][2-body]") {

  std::vector<Shell> obs{
      // pseudorandom d
      Shell{{1.0, 3.0}, {{2, true, {1.0, 0.3}}}, {{0.0, 0.0, 0.0}}},
      // pseudorandom d
      Shell{{2.0, 5.0}, {{2, true, {1.0, 0.2}}}, {{1.0, 1.0, 1.0}}},
//    tight functions are problematic with Gaussian fits. The errors in some integrals involving the functions below can be huge
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
//      Shell{{11720.0000000, 1759.0000000, 400.8000000, 113.7000000, 37.0300000,
//             13.2700000, 5.0250000, 1.0130000},
//            {{0,
//              false,
//              {0.0007100, 0.0054700, 0.0278370, 0.1048000, 0.2830620, 0.4487190,
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
//  std::vector<std::pair<double, double>> cgtg_params{{0.2209,0.3144},{1.004,0.3037},{3.622,0.1681},{12.16,0.09811},{45.87,0.06024},{254.4,0.03726}};
  // "precise" fit of exp(-r12) on [0,10)
  std::vector<std::pair<double,double>> cgtg_params =
    {{0.10535330565471572,0.08616353459042002},{0.22084823136587992,0.08653979627551414},{0.3543431104992702,0.08803697599356214},{0.48305514665749105,0.09192519612306953},{0.6550035700167584,0.10079776426873248},{1.1960917050801643,0.11666110644901695},{2.269278814810891,0.14081371547404428},{5.953990617813977,0.13410255216448014},{18.31911063199608,0.0772095196191394},{66.98443868169818,0.049343985939540556},{367.24137290439205,0.03090625839896873},{5.655142311118115,-0.017659052507938647}};

  if (LIBINT2_MAX_AM_eri >= max_l) {
    for (int k = -1; k <= 0; ++k) {
      auto engine = Engine(k == 0 ? Operator::stg : Operator::yukawa, max_nprim, max_l);
      const auto scale = 2.3;
      engine.set_params(1.0).prescale_by(scale);
      REQUIRE(engine.prescaled_by() == scale);
      const auto &results = engine.results();
      auto engine_ref =
          Engine(k == 0 ? Operator::cgtg : Operator::cgtg_x_coulomb, max_nprim, max_l);
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
                  REQUIRE(results[0][i]/scale ==
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
TEST_CASE_METHOD(libint2::unit::DefaultFixture, "bra-ket permutation", "[engine][2-body]") {
  using libint2::Shell;
  using libint2::BasisSet;
  using libint2::Operator;
  using libint2::Engine;
  std::vector<Shell> shells;
  shells.push_back({
      {0.8378385011e+02, 0.1946956493e+02, 0.6332106784e+01}, // exponents
      {// P shell, spherical=false, contraction coefficients
       {1, false, {0.1559162750, 0.6076837186, 0.3919573931}}},
      {{0, 0, 0}} // origin coordinates
  });
  shells.push_back({
      {0.2964817927e+01, 0.9043639676e+00, 0.3489317337e+00}, // exponents
      {// D shell, spherical=false, contraction coefficients
       {2, false, {0.2197679508, 0.6555473627, 0.2865732590}}},
      {{0, 0, 0}} // origin coordinates
  });

  auto obs = BasisSet();
  for (auto&& shell : shells) {
    obs.push_back(std::move(shell));
  }
  auto shell2bf = BasisSet::compute_shell2bf(obs);

  auto engine =
      Engine(libint2::Operator::coulomb, libint2::max_nprim(obs),
                      libint2::max_l(obs), 0);
  auto engine_kb =
      Engine(libint2::Operator::coulomb, libint2::max_nprim(obs),
             libint2::max_l(obs), 0);
  // Force uniform normalised Cartesian functions
  engine.set(libint2::CartesianShellNormalization::uniform);
  engine_kb.set(libint2::CartesianShellNormalization::uniform);
  const auto &buf = engine.results();
  const auto &buf_kb = engine_kb.results();

  for (auto s1 = 0; s1 != obs.size(); ++s1) {
    const auto bf1_first = shell2bf[s1]; // first basis function in this shell
    auto n1 = obs[s1].size();            // number of basis functions in shell 1
    for (auto s2 = 0; s2 != obs.size(); ++s2) {
      const auto bf2_first = shell2bf[s2]; // first basis function in this shell
      auto n2 = obs[s2].size(); // number of basis functions in shell 2
      for (auto s3 = 0; s3 != obs.size(); ++s3) {
        const auto bf3_first =
            shell2bf[s3];         // first basis function in this shell
        auto n3 = obs[s3].size(); // number of basis functions in shell 3
        for (auto s4 = 0; s4 != obs.size(); ++s4) {
          const auto bf4_first =
              shell2bf[s4];         // first basis function in this shell
          auto n4 = obs[s4].size(); // number of basis functions in shell 4
          engine.compute2<Operator::coulomb, libint2::BraKet::xx_xx, 0>(
              obs[s1], obs[s2], obs[s3], obs[s4]);
          const auto *buf_1234 = buf[0];
          engine_kb.compute2<Operator::coulomb, libint2::BraKet::xx_xx, 0>(
              obs[s3], obs[s4], obs[s1], obs[s2]);
          const auto *buf_3412 = buf_kb[0];
          for (auto f1 = 0ul, f1234 = 0ul; f1 != n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for (auto f2 = 0ul; f2 != n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for (auto f3 = 0ul; f3 != n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for (auto f4 = 0ul; f4 != n4; ++f4, ++f1234) {
                  const auto bf4 = f4 + bf4_first;
                  const auto integral = buf_1234[f1234];
                  const auto f3412 = ((f3 * n4 + f4)*n1 + f1)*n2 + f2;
                  const auto integral_kb = buf_3412[f3412];
                  REQUIRE(std::abs(integral-integral_kb) < 1e-14);
                }
              }
            }
          }
        }
      }
    }
  }
}

