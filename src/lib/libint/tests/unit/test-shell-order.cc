#include "catch.hpp"
#include "fixture.h"

#include <libint2/cgshell_ordering.h>
#include <libint2/shgshell_ordering.h>


TEST_CASE_METHOD(libint2::unit::DefaultFixture, "solid harmonics", "[engine]") {

  std::vector<Shell> obs{
      Shell{{1.0, 3.0}, {{1, true, {1.0, 0.3}}}, {{0.0, 0.0, 0.0}}},
      Shell{{1.0, 3.0},
            {{1, false, {1.0, 0.3}}},
            {{0.0, 0.0, 0.0}}}, // same as 0, but Cartesian
      Shell{{2.0, 5.0}, {{2, true, {1.0, 0.2}}}, {{1.0, 1.0, 1.0}}}};

  // test that integrals with p solid harmonics are correctly computed
  // N.B. since by default BasisSet makes p's cartesian this is rarely (if at all) tested in production
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
      FOR_CART(ix, iy, iz, 1)
      int lm;
      if (ix == 1) {
        lm = INT_SOLIDHARMINDEX(1, 1);
      } else if (iy == 1) {
        lm = INT_SOLIDHARMINDEX(1, -1);
      } else if (iz == 1) {
        lm = INT_SOLIDHARMINDEX(1, 0);
      }

      for (int f2 = 0; f2 != 5; ++f2) {
        auto cart_f1f2 = xyz * 5 + f2;
        auto sph_f1f2 = lm * 5 + f2;
        REQUIRE(PD_ints[sph_f1f2] == Approx(pD_ints[cart_f1f2]));
      }

      ++xyz;
      END_FOR_CART
    }

#endif // LIBINT2_SUPPORT_ONEBODY
  }  // SECTION("solid harmonic p")

}
