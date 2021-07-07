//
// Created by Eduard Valeyev on 7/28/20.
//

#include <array>
#include <iostream>

#include "catch.hpp"
#include <libint2/config.h>
#include <libint2/numeric.h>
#include <libint2/engine.h>
#include <test_eri/eri.h>

#ifdef LIBINT_HAS_MPFR
TEST_CASE("ERI reference values", "[2-body][precision]") {

  {
    /// Xiaosong's test case:
    /// \f[
    /// \sum_i d/dAi d/dBi (D_xz D_xz|D_xz D_xz)
    /// \f]
    /// where \f$ D_xy \f$ is primitive D function with exponent 15000

    const LIBINT2_REF_REALTYPE alpha = 15000;
    const LIBINT2_REF_REALTYPE O[] = {0,0,0};
    LIBINT2_REF_REALTYPE ref_value = 0;
    for (int xyz = 0; xyz != 3; ++xyz) {
      std::array<unsigned int, 12> deriv_idx = {
          {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
      deriv_idx[xyz] = 1;
      deriv_idx[3 + xyz] = 1;
      const int norm_flag = 1;
      ref_value += eri(deriv_idx.data(), 1, 1, 0, alpha, O, 1, 1, 0, alpha, O, 1, 1, 0, alpha, O, 1, 1, 0, alpha, O,
                   norm_flag);
    }
//    std::cout << std::setprecision(20) << "ref_value=" << ref_value << std::endl;
//    mp_exp_t exp;
//    std::cout << std::setprecision(20) << "ref_value=" << ref_value.get_str(exp) << "e" << exp << std::endl;

    using namespace libint2;
    Shell d{{15000}, {{2, true, {1.0}}}, {{0.0, 0.0, 0.0}}};
    Engine engine(Operator::coulomb, 1, 2, 2);
    const auto &results = engine.results();
    engine.compute(d,d,d,d);
    /*
     * Use this Mathematica:
     * tri[i_, j_] := -12 - 1/2 (-24 + Min[i, j]) (1 + Min[i, j]) + Max[i, j];
     * tri[0, 3]
     * tri[1, 4]
     * tri[2, 5] 
     */
    auto value = results[3][0]+results[15][0]+results[26][0];
//    std::cout << "value=" << value << std::endl;

    using std::abs;
    CHECK(abs((ref_value - value)/ref_value) < 2e-15);
  }
}
#endif  // LIBINT_HAS_MPFR

