#include "catch.hpp"
#include "fixture.h"

TEST_CASE_METHOD(libint2::unit::DefaultFixture, "electrostatic potential", "[engine][1-body]") {
#if defined(LIBINT2_SUPPORT_ONEBODY)
  if (LIBINT_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD)
    return;

  std::vector<Shell> obs{
      Shell{{1.0, 3.0}, {{2, true, {1.0, 0.3}}}, {{0.0, 0.0, 0.0}}},
      Shell{{2.0, 5.0}, {{2, true, {1.0, 0.2}}}, {{1.0, 1.0, 1.0}}}};
  {
    const auto lmax = std::min(3, LIBINT2_MAX_AM_elecpot);
    if (lmax >= 2) {
      auto engine = Engine(Operator::nuclear, 2, lmax);
      engine.set_params(make_point_charges(atoms));

      engine.compute(obs[0], obs[0]);
      {
        std::vector<double> shellset_ref = {
            -1.238239259091998e+01, 0.000000000000000e+00,
            0.000000000000000e+00,  -5.775996163160049e-02,
            0.000000000000000e+00,  0.000000000000000e+00,
            -1.301230978657952e+01, -6.796143730068988e-02,
            0.000000000000000e+00,  1.139389632827834e-01,
            0.000000000000000e+00,  -6.796143730068988e-02,
            -1.343732979153083e+01, 0.000000000000000e+00,
            -1.478824785355970e-02, -5.775996163160049e-02,
            0.000000000000000e+00,  0.000000000000000e+00,
            -1.284475452992947e+01, 0.000000000000000e+00,
            0.000000000000000e+00,  1.139389632827834e-01,
            -1.478824785355970e-02, 0.000000000000000e+00,
            -1.241040347301479e+01};
        for (int i = 0; i != 25; ++i) {
          REQUIRE(engine.results()[0][i] == Approx(shellset_ref[i]));
        }
      }

      engine.compute(obs[0], obs[1]);
      {
        std::vector<double> shellset_ref = {
            -4.769186621041819e-01, -9.303619356400431e-01,
            -1.559058302243514e+00, -9.290824121864600e-01,
            -5.835786921473129e-04, -1.159266418436018e+00,
            -3.770080831197964e-01, 9.572841308198474e-01,
            -8.291498398421207e-01, -1.663667687168316e+00,
            -2.171951144148577e+00, 1.074249956874296e+00,
            2.128355904665372e+00,  1.074590109905394e+00,
            -3.485163651594458e-03, -1.160865205880651e+00,
            -8.344173649626901e-01, 9.566621490332916e-01,
            -3.760919234260182e-01, 1.660514988916377e+00,
            -1.120272634615116e-03, -1.385603731947886e+00,
            -2.105750177166632e-03, 1.380654897976564e+00,
            2.115041199099945e+00};
        for (int i = 0; i != 25; ++i) {
          REQUIRE(engine.results()[0][i] == Approx(shellset_ref[i]));
        }
      }
    }
  }

#endif // LIBINT2_SUPPORT_ONEBODY
}
