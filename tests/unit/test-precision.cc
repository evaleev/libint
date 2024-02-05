/*
 *  Copyright (C) 2004-2024 Edward F. Valeev
 *
 *  This file is part of Libint library.
 *
 *  Libint library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <libint2/config.h>
#include <libint2/engine.h>
#include <libint2/numeric.h>

#include <array>
#include <iostream>

#include "catch.hpp"
#if defined(NO_LIBINT_COMPILER_CODE)
#include "../eri/eri.h"
#else
#include <test_eri/eri.h>
#endif

#ifdef LIBINT_HAS_MPFR
TEST_CASE("ERI reference values", "[2-body][precision]") {
  {
    /// Xiaosong's test case:
    /// \f[
    /// \sum_i d/dAi d/dBi (D_xz D_xz|D_xz D_xz)
    /// \f]
    /// where \f$ D_xy \f$ is primitive D function with exponent 15000

    const LIBINT2_REF_REALTYPE alpha = 15000;
    const LIBINT2_REF_REALTYPE O[] = {0, 0, 0};
    LIBINT2_REF_REALTYPE ref_value = 0;
    for (int xyz = 0; xyz != 3; ++xyz) {
      std::array<unsigned int, 12> deriv_idx = {
          {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
      deriv_idx[xyz] = 1;
      deriv_idx[3 + xyz] = 1;
      const int norm_flag = 1;
      ref_value += eri(deriv_idx.data(), 1, 1, 0, alpha, O, 1, 1, 0, alpha, O,
                       1, 1, 0, alpha, O, 1, 1, 0, alpha, O, norm_flag);
    }
    //    std::cout << std::setprecision(20) << "ref_value=" << ref_value <<
    //    std::endl; mp_exp_t exp; std::cout << std::setprecision(20) <<
    //    "ref_value=" << ref_value.get_str(exp) << "e" << exp << std::endl;

    using namespace libint2;
    Shell d{{15000}, {{2, true, {1.0}}}, {{0.0, 0.0, 0.0}}};
    Engine engine(Operator::coulomb, 1, 2, 2);
    const auto& results = engine.results();
    engine.compute(d, d, d, d);
    /*
     * Use this Mathematica:
     * tri[i_, j_] := -12 - 1/2 (-24 + Min[i, j]) (1 + Min[i, j]) + Max[i, j];
     * tri[0, 3]
     * tri[1, 4]
     * tri[2, 5]
     */
    auto value = results[3][0] + results[15][0] + results[26][0];
    //    std::cout << "value=" << value << std::endl;

    using std::abs;
    CHECK(abs((ref_value - value) / ref_value) < 2e-15);
  }
}
#endif  // LIBINT_HAS_MPFR

#if defined(LIBINT2_SUPPORT_ERI)
TEST_CASE("2-body integrals precision", "[2-body][precision]") {
  using namespace libint2;

  const auto original_default_screening_method =
      libint2::default_screening_method();
  libint2::default_screening_method(libint2::ScreeningMethod::Conservative);

  const auto ry = 4.0;
  const auto rz = 20.0;
  const std::string obs_name = "sto-3g";
  const std::string dfbs_name = "cc-pvdz";
  //  const auto ry = 5.0;
  //  const auto rz = 40.0;
  //  const std::string obs_name = "cc-pvtz";
  //  const std::string dfbs_name = "cc-pvqz-ri";
  std::stringstream sstr;
  sstr << "3\n\nNe 0 0 0\nNe 0 0 " << rz << "\nNe 0 " << ry << " 0\n";
  auto atoms = libint2::read_dotxyz(sstr);
  BasisSet obs(obs_name, atoms);
  BasisSet dfbs(dfbs_name, atoms);
  const auto max_nprim = std::max(obs.max_nprim(), dfbs.max_nprim());
  const auto max_l = std::max(obs.max_l(), dfbs.max_l());

  for (auto eps : {1e-8, 1e-10, 1e-12, 1e-14, 1e-16}) {
    for (auto oper : {Operator::coulomb}) {
      for (auto c : {2, 3, 4}) {
#if defined(LIBINT2_SUPPORT_ERI2)
        if (max_l > LIBINT2_MAX_AM_2eri) continue;
#else
        if (c == 2) continue;
#endif
#if defined(LIBINT2_SUPPORT_ERI3)
        if (max_l > LIBINT2_MAX_AM_3eri) continue;
#else
        if (c == 3) continue;
#endif
#if defined(LIBINT2_SUPPORT_ERI)
        if (max_l > LIBINT2_MAX_AM_eri) continue;
#else
        if (c == 4) continue;
#endif

        auto bases = (c == 4)
                         ? std::vector<BasisSet>{obs, obs, obs, obs}
                         : ((c == 3) ? std::vector<BasisSet>{dfbs, obs, obs}
                                     : std::vector<BasisSet>{dfbs, dfbs});

        auto braket = (c == 4) ? BraKet::xx_xx
                               : ((c == 3) ? BraKet::xs_xx : BraKet::xs_xs);

        Engine ref_engine(oper, max_nprim, max_l, 0);
        ref_engine.set_precision(0.0);
        ref_engine.set(braket);
        Engine eps_engine(oper, max_nprim, max_l, 0);
        eps_engine.set_precision(eps);
        eps_engine.set(braket);

        auto test = [&ref_engine, &eps_engine](
                        const std::vector<const Shell*>& shells,
                        double precision) {
          const auto c = shells.size();
          assert(c == 4 || c == 3 || c == 2);

          const auto ref_ints =
              (c == 4)
                  ? ref_engine.compute(*shells[0], *shells[1], *shells[2],
                                       *shells[3])
                  : ((c == 3) ? ref_engine.compute(*shells[0], *shells[1],
                                                   *shells[2])
                              : ref_engine.compute(*shells[0], *shells[1]));
          assert(ref_ints[0] != nullptr);

          const auto eps_ints =
              (c == 4)
                  ? eps_engine.compute(*shells[0], *shells[1], *shells[2],
                                       *shells[3])
                  : ((c == 3) ? eps_engine.compute(*shells[0], *shells[1],
                                                   *shells[2])
                              : eps_engine.compute(*shells[0], *shells[1]));

          const auto nf = std::accumulate(shells.begin(), shells.end(), 1,
                                          [](int nf, const Shell* shell_ptr) {
                                            return nf * shell_ptr->size();
                                          });

          auto max_engine_screening_error = 0.;
          LIBINT_MAYBE_UNUSED bool engine_precision_too_low = false;
          const auto max_allowed_exceening_error_factor = 2;
          for (auto f = 0; f != nf; ++f) {
            auto& ref_v = ref_ints[0][f];
            const auto engine_screening_error =
                (eps_ints[0] != nullptr) ? std::abs(eps_ints[0][f] - ref_v)
                                         : std::abs(ref_v);
            CHECK(engine_screening_error <=
                  max_allowed_exceening_error_factor * precision);
            if (engine_screening_error > precision) {
              max_engine_screening_error =
                  std::max(engine_screening_error, max_engine_screening_error);
            }
            if (engine_screening_error >
                max_allowed_exceening_error_factor * precision) {
              engine_precision_too_low = true;
            }
          }
        };

        std::vector<const Shell*> shells(c);
        for (auto&& sh0 : bases[0]) {
          shells[0] = &sh0;
          for (auto&& sh1 : bases[0]) {
            shells[1] = &sh1;
            if (c > 2) {
              for (auto&& sh2 : bases[2]) {
                shells[2] = &sh2;
                if (c > 3) {
                  for (auto&& sh3 : bases[3]) {
                    shells[3] = &sh3;
                    test(shells, eps);
                  }
                } else {
                  test(shells, eps);
                }
              }
            } else {
              test(shells, eps);
            }
          }
        }
      }
    }
  }

  libint2::default_screening_method(original_default_screening_method);
}
#endif  // defined(LIBINT2_SUPPORT_ERI)
