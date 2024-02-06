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

#include "catch.hpp"
#include "fixture.h"

std::tuple<int, int, int, int> permute4_1234(int i1, int i2, int i3, int i4) {
  return std::make_tuple(i1, i2, i3, i4);
};
std::tuple<int, int, int, int> permute4_2134(int i1, int i2, int i3, int i4) {
  return std::make_tuple(i2, i1, i3, i4);
};
std::tuple<int, int, int, int> permute4_1243(int i1, int i2, int i3, int i4) {
  return std::make_tuple(i1, i2, i4, i3);
};
std::tuple<int, int, int, int> permute4_2143(int i1, int i2, int i3, int i4) {
  return std::make_tuple(i2, i1, i4, i3);
};
std::tuple<int, int, int, int> permute4_3412(int i1, int i2, int i3, int i4) {
  return std::make_tuple(i3, i4, i1, i2);
};
std::tuple<int, int, int, int> permute4_3421(int i1, int i2, int i3, int i4) {
  return std::make_tuple(i3, i4, i2, i1);
};
std::tuple<int, int, int, int> permute4_4312(int i1, int i2, int i3, int i4) {
  return std::make_tuple(i4, i3, i1, i2);
};
std::tuple<int, int, int, int> permute4_4321(int i1, int i2, int i3, int i4) {
  return std::make_tuple(i4, i3, i2, i1);
};

std::tuple<int, int, int> permute3_123(int i1, int i2, int i3) {
  return std::make_tuple(i1, i2, i3);
};
std::tuple<int, int, int> permute3_213(int i1, int i2, int i3) {
  return std::make_tuple(i2, i1, i3);
};
std::tuple<int, int, int> permute3_132(int i1, int i2, int i3) {
  return std::make_tuple(i1, i3, i2);
};
std::tuple<int, int, int> permute3_231(int i1, int i2, int i3) {
  return std::make_tuple(i2, i3, i1);
};
std::tuple<int, int, int> permute3_312(int i1, int i2, int i3) {
  return std::make_tuple(i3, i1, i2);
};
std::tuple<int, int, int> permute3_321(int i1, int i2, int i3) {
  return std::make_tuple(i3, i2, i1);
};

std::tuple<int, int> permute2_12(int i1, int i2) {
  return std::make_tuple(i1, i2);
};
std::tuple<int, int> permute2_21(int i1, int i2) {
  return std::make_tuple(i2, i1);
};

/// takes a composite (upper-triangle) 2nd deriv index and returns
/// {smaller,greater} of the 2 derivative indices
std::tuple<int, int> split_deriv2(int natoms, int d) {
  const int n = natoms * 3;
  const int dmin =
      std::ceil((-1 + 2 * n - std::sqrt(-7. - 8 * d + 4 * n + 4 * n * n)) / 2);
  const int dmax = d - (dmin * ((2 * n) - dmin - 1)) / 2;
  return std::make_tuple(dmin, dmax);
}

int merge_deriv2(int natoms, int d1, int d2) {
  const int n = natoms * 3;  // # of order-1 derivatives
  const int dmin = std::min(d1, d2);
  const int dmax = std::max(d1, d2);
  return (dmin * ((2 * n) - dmin - 1)) / 2 + dmax;
}

template <unsigned int deriv_order>
void validate4(const BasisSet& obs, const std::vector<Atom>& atoms) {
#ifdef LIBINT2_SUPPORT_ERI
  constexpr int ncenters = 4;
  const auto max_l = obs.max_l();
  if (deriv_order > LIBINT2_DERIV_ERI_ORDER) return;
  switch (deriv_order) {
    case 0:
      if (max_l > LIBINT2_MAX_AM_eri) return;
      break;
    case 1:
#if LIBINT2_DERIV_ERI_ORDER > 0
      if (max_l > LIBINT2_MAX_AM_eri1) return;
#endif
      break;
    case 2:
#if LIBINT2_DERIV_ERI_ORDER > 1
      if (max_l > LIBINT2_MAX_AM_eri2) return;
#endif
      break;
    default:
      abort();
  }

  const auto abs_precision = deriv_order == 0 ? 1e-13 : 1e-12;

  Engine engine_ref(Operator::coulomb, obs.max_nprim(), obs.max_l(),
                    deriv_order);
  Engine engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), deriv_order);

  const auto nderiv_shellset = libint2::num_geometrical_derivatives(
      4, deriv_order);  // # of derivs for each shell quartet

  auto shell2bf = obs.shell2bf();  // maps shell index to basis function index
  auto shell2atom = obs.shell2atom(atoms);

  for (auto s1 = 0; s1 != obs.size(); ++s1) {
    for (auto s2 = 0; s2 <= s1; ++s2) {
      for (auto s3 = 0; s3 <= s1; ++s3) {
        auto s4_max = (s1 == s3) ? s2 : s3;
        for (auto s4 = 0; s4 <= s4_max; ++s4) {
          engine_ref.compute(obs[s1], obs[s2], obs[s3], obs[s4]);
          const auto& shellset_ref =
              engine_ref.results();  // location of the computed integrals
          if (shellset_ref[0] == nullptr)
            continue;  // nullptr returned if the entire shell-set was screened
                       // out

          auto n1 = obs[s1].size();  // number of basis functions in first shell
          auto n2 =
              obs[s2].size();  // number of basis functions in second shell
          auto n3 = obs[s3].size();  // number of basis functions in third shell
          auto n4 =
              obs[s4].size();  // number of basis functions in fourth shell

          auto validate = [&](std::tuple<int, int, int, int> (&permuter)(
                              int, int, int, int)) {
            int ss1, ss2, ss3, ss4;
            std::tie(ss1, ss2, ss3, ss4) = permuter(s1, s2, s3, s4);

            engine.compute(obs[ss1], obs[ss2], obs[ss3], obs[ss4]);
            const auto& shellset =
                engine.results();  // location of the computed integrals

            REQUIRE(shellset[0] != nullptr);
            for (int d = 0; d != nderiv_shellset; ++d) {
              auto compare = [&](int dd) {
                int nn1, nn2, nn3, nn4;
                std::tie(nn1, nn2, nn3, nn4) = permuter(n1, n2, n3, n4);
                for (auto f1 = 0; f1 != n1; ++f1) {
                  for (auto f2 = 0; f2 != n2; ++f2) {
                    for (auto f3 = 0; f3 != n3; ++f3) {
                      for (auto f4 = 0; f4 != n4; ++f4) {
                        int ff1, ff2, ff3, ff4;
                        std::tie(ff1, ff2, ff3, ff4) = permuter(f1, f2, f3, f4);

                        auto f1234 = ((f1 * n2 + f2) * n3 + f3) * n4 + f4;
                        auto ff1234 =
                            ((ff1 * nn2 + ff2) * nn3 + ff3) * nn4 + ff4;

                        REQUIRE(shellset[dd][ff1234] ==
                                Approx(shellset_ref[d][f1234])
                                    .margin(abs_precision));
                      }
                    }
                  }
                }
              };

              auto permute_deriv = [&](const int d) {
                assert(ncenters == 4);
                int dd = -1;
                const auto xyz = d % 3;
                const auto atom = d / 3;
                int a[4] = {0, 0, 0, 0};
                a[atom] = 1;
                int aa[4];
                std::tie(aa[0], aa[1], aa[2], aa[3]) =
                    permuter(a[0], a[1], a[2], a[3]);
                if (aa[0])
                  dd = xyz;
                else if (aa[1])
                  dd = 3 + xyz;
                else if (aa[2])
                  dd = 6 + xyz;
                else if (aa[3])
                  dd = 9 + xyz;
                else
                  assert(false);
                return dd;
              };

              // permuting geometric deriv index is order-specific
              int dd;
              switch (deriv_order) {
                case 0:
                  dd = 0;
                  break;

                case 1: {
                  dd = permute_deriv(d);
                } break;

                case 2: {
                  int d1, d2;
                  std::tie(d1, d2) = split_deriv2(ncenters, d);
                  auto dd1 = permute_deriv(d1);
                  auto dd2 = permute_deriv(d2);
                  dd = merge_deriv2(ncenters, dd1, dd2);
                } break;

                default:
                  assert(false);  // can't compare yet deriv_order > 2
              }
              compare(dd);

            }  // d
          };

          validate(permute4_2134);
          validate(permute4_1243);
          validate(permute4_2143);
          validate(permute4_3412);
          validate(permute4_3421);
          validate(permute4_4312);
          validate(permute4_4321);
        }
      }
    }
  }
#endif  // LIBINT2_SUPPORT_ERI
}

template <unsigned int deriv_order>
void validate3(const BasisSet& obs, const BasisSet& dfbs,
               const std::vector<Atom>& atoms) {
#if defined(LIBINT2_SUPPORT_ERI) && defined(LIBINT2_SUPPORT_ERI3)
  constexpr int ncenters = 3;
  const auto max_l = std::max(obs.max_l(), dfbs.max_l());
  if (deriv_order > LIBINT2_DERIV_ERI_ORDER ||
      deriv_order > LIBINT2_DERIV_ERI3_ORDER)
    return;
  const auto xsxx = LIBINT_SHELL_SET == LIBINT_SHELL_SET_STANDARD;
  if (!xsxx) return;  // not yet implemented
  const auto abs_precision = deriv_order == 1 ? 1e-12 : 1e-13;

  const auto max_nprim = std::max(obs.max_nprim(), dfbs.max_nprim());
  Engine engine_ref(Operator::coulomb, max_nprim, max_l, deriv_order);
  engine_ref.set(BraKet::xx_xx);
  Engine engine(Operator::coulomb, max_nprim, max_l, deriv_order);
  engine.set(xsxx ? BraKet::xs_xx : BraKet::xx_xs);

  const auto nderiv_shellset_ref = libint2::num_geometrical_derivatives(
      4, deriv_order);  // # of derivs for reference shell *quartet*

  auto shell2bf = obs.shell2bf();  // maps shell index to basis function index
  auto shell2atom = obs.shell2atom(atoms);
  auto shell2bf_df =
      dfbs.shell2bf();  // maps shell index to basis function index
  auto shell2atom_df = dfbs.shell2atom(atoms);

  for (auto s1 = 0; s1 != dfbs.size(); ++s1) {
    for (auto s2 = 0; s2 != obs.size(); ++s2) {
      for (auto s3 = 0; s3 != obs.size(); ++s3) {
        assert(xsxx);

        // skip if angular momenta are too high
        const auto max_orb_l = std::max(obs[s2].contr[0].l, obs[s3].contr[0].l);
        const auto max_l = std::max(dfbs[s1].contr[0].l, max_orb_l);
        auto max_l_exceeded = false;
        switch (deriv_order) {
          case 0:
            if (max_l > LIBINT2_MAX_AM_eri || max_l > LIBINT2_MAX_AM_3eri)
              max_l_exceeded = true;
#ifdef LIBINT2_CENTER_DEPENDENT_MAX_AM_3eri
            if (max_orb_l > LIBINT2_MAX_AM_default) max_l_exceeded = true;
#endif
            break;

          case 1:
#if LIBINT2_DERIV_ERI_ORDER > 0 && LIBINT2_DERIV_ERI3_ORDER > 0
            if (max_l > LIBINT2_MAX_AM_eri1 || max_l > LIBINT2_MAX_AM_3eri1)
              max_l_exceeded = true;
#ifdef LIBINT2_CENTER_DEPENDENT_MAX_AM_3eri1
            if (max_orb_l > LIBINT2_MAX_AM_default1) max_l_exceeded = true;
#endif
#endif
            break;
          case 2:
#if LIBINT2_DERIV_ERI_ORDER > 1 && LIBINT2_DERIV_ERI3_ORDER > 1
            if (max_l > LIBINT2_MAX_AM_eri2 || max_l > LIBINT2_MAX_AM_3eri2)
              max_l_exceeded = true;
#ifdef LIBINT2_CENTER_DEPENDENT_MAX_AM_3eri1
            if (max_orb_l > LIBINT2_MAX_AM_default2) max_l_exceeded = true;
#endif
#endif
            break;
          default:
            abort();
        }
        if (max_l_exceeded) continue;

        engine_ref.compute(dfbs[s1], Shell::unit(), obs[s2], obs[s3]);

        const auto& shellset_ref =
            engine_ref.results();  // location of the computed integrals
        if (shellset_ref[0] == nullptr)
          continue;  // nullptr returned if the entire shell-set was screened
                     // out

        auto n1 = dfbs[s1].size();  // number of basis functions in first shell
        auto n2 = obs[s2].size();   // number of basis functions in second shell
        auto n3 = obs[s3].size();   // number of basis functions in third shell

        auto validate = [&](std::tuple<int, int, int> (&permuter)(int, int,
                                                                  int)) {
          int ss1, ss2, ss3;
          std::tie(ss1, ss2, ss3) = permuter(s1, s2, s3);

          assert(xsxx);
          engine.compute(dfbs[ss1], obs[ss2], obs[ss3]);
          const auto& shellset =
              engine.results();  // location of the computed integrals
          REQUIRE(shellset[0] != nullptr);

          for (int d = 0; d != nderiv_shellset_ref; ++d) {
            auto compare = [&](int dd) {
              assert(xsxx);  // xxxs is not implemented
              int nn1, nn2, nn3;
              std::tie(nn1, nn2, nn3) = permuter(n1, n2, n3);
              for (auto f1 = 0; f1 != n1; ++f1) {
                for (auto f2 = 0; f2 != n2; ++f2) {
                  for (auto f3 = 0; f3 != n3; ++f3) {
                    int ff1, ff2, ff3;
                    std::tie(ff1, ff2, ff3) = permuter(f1, f2, f3);

                    auto f123 = (f1 * n2 + f2) * n3 + f3;
                    auto ff123 = (ff1 * nn2 + ff2) * nn3 + ff3;

                    REQUIRE(
                        shellset[dd][ff123] ==
                        Approx(shellset_ref[d][f123]).margin(abs_precision));
                  }
                }
              }
            };

            auto permute_deriv = [&](const int d) {
              assert(ncenters == 3);
              int dd;
              const auto xyz = d % 3;
              const auto atom = d / 3;

              // skip derivatives w.r.t. the dummy atom
              assert(xsxx);
              if (atom == 1) return -1;

              int a[3] = {0, 0, 0};
              a[atom == 0 ? 0 : atom - 1] = 1;
              int aa[3];
              assert(xsxx);
              std::tie(aa[0], aa[1], aa[2]) = permuter(a[0], a[1], a[2]);
              if (aa[0])
                dd = xyz;
              else if (aa[1])
                dd = 3 + xyz;
              else if (aa[2])
                dd = 6 + xyz;
              else
                assert(false);
              return dd;
            };

            // permuting geometric deriv index is order-specific
            int dd;
            switch (deriv_order) {
              case 0:
                dd = 0;
                break;

              case 1: {
                dd = permute_deriv(d);
                if (dd == -1) return;  // skip unneeded derivatives
              } break;

              case 2: {
                int d1, d2;
                std::tie(d1, d2) = split_deriv2(ncenters, d);
                auto dd1 = permute_deriv(d1);
                auto dd2 = permute_deriv(d2);
                if (dd1 == -1 || dd2 == -1)
                  return;  // skip unneeded derivatives
                dd = merge_deriv2(ncenters, dd1, dd2);
              } break;

              default:
                assert(false);  // can't compare yet deriv_order > 2
            }
            compare(dd);

          }  // d
        };

        validate(permute3_123);
        validate(permute3_132);
      }
    }
  }
#endif  // LIBINT2_SUPPORT_ERI3
}

template <unsigned int deriv_order>
void validate2(const BasisSet& obs, const std::vector<Atom>& atoms) {
#if defined(LIBINT2_SUPPORT_ERI) && defined(LIBINT2_SUPPORT_ERI2)
  constexpr int ncenters = 2;
  if (deriv_order > LIBINT2_DERIV_ERI_ORDER ||
      deriv_order > LIBINT2_DERIV_ERI2_ORDER)
    return;
  switch (deriv_order) {
    case 0:
      if (obs.max_l() > LIBINT2_MAX_AM_eri || obs.max_l() > LIBINT2_MAX_AM_2eri)
        return;
      break;
    case 1:
#if LIBINT2_DERIV_ERI_ORDER > 0 && LIBINT2_DERIV_ERI2_ORDER > 0
      if (obs.max_l() > LIBINT2_MAX_AM_eri1 ||
          obs.max_l() > LIBINT2_MAX_AM_2eri1)
        return;
#endif
      break;
    case 2:
#if LIBINT2_DERIV_ERI_ORDER > 1 && LIBINT2_DERIV_ERI2_ORDER > 1
      if (obs.max_l() > LIBINT2_MAX_AM_eri2 ||
          obs.max_l() > LIBINT2_MAX_AM_2eri2)
        return;
#endif
      break;
    default:
      abort();
  }
  const auto xsxs = LIBINT_SHELL_SET == LIBINT_SHELL_SET_STANDARD;
  if (!xsxs) return;  // not yet implemented
  const auto abs_precision = deriv_order == 1 ? 1e-12 : 1e-13;

  Engine engine_ref(Operator::coulomb, obs.max_nprim(), obs.max_l(),
                    deriv_order);
  engine_ref.set(BraKet::xx_xx);
  Engine engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), deriv_order);
  assert(xsxs);
  engine.set(BraKet::xs_xs);

  const auto nderiv_shellset_ref = libint2::num_geometrical_derivatives(
      4, deriv_order);  // # of derivs for reference shell *quartet*

  auto shell2bf = obs.shell2bf();  // maps shell index to basis function index
  auto shell2atom = obs.shell2atom(atoms);

  for (auto s1 = 0; s1 != obs.size(); ++s1) {
    for (auto s2 = 0; s2 != obs.size(); ++s2) {
      assert(xsxs);
      engine_ref.compute(obs[s1], Shell::unit(), obs[s2], Shell::unit());

      const auto& shellset_ref =
          engine_ref.results();  // location of the computed integrals
      if (shellset_ref[0] == nullptr)
        continue;  // nullptr returned if the entire shell-set was screened out

      auto n1 = obs[s1].size();  // number of basis functions in first shell
      auto n2 = obs[s2].size();  // number of basis functions in second shell

      auto validate = [&](std::tuple<int, int> (&permuter)(int, int)) {
        int ss1, ss2;
        std::tie(ss1, ss2) = permuter(s1, s2);

        assert(xsxs);
        engine.compute(obs[ss1], obs[ss2]);
        const auto& shellset =
            engine.results();  // location of the computed integrals

        REQUIRE(shellset[0] != nullptr);
        for (int d = 0; d != nderiv_shellset_ref; ++d) {
          auto compare = [&](int dd) {
            assert(xsxs);  // xxxs is not implemented
            int nn1, nn2;
            std::tie(nn1, nn2) = permuter(n1, n2);
            for (auto f1 = 0; f1 != n1; ++f1) {
              for (auto f2 = 0; f2 != n2; ++f2) {
                int ff1, ff2;
                std::tie(ff1, ff2) = permuter(f1, f2);

                auto f12 = f1 * n2 + f2;
                auto ff12 = ff1 * nn2 + ff2;

                if (std::abs(shellset[dd][ff12] - shellset_ref[d][f12]) >=
                    abs_precision) {
                  const auto value = shellset[dd][ff12];
                  const auto value_ref = shellset_ref[d][f12];
                  const auto value_minus_value_ref = value - value_ref;
                  std::cout << "|value - value_ref| = " << std::scientific
                            << value_minus_value_ref << " value = " << value
                            << std::endl;
                }
                REQUIRE(shellset[dd][ff12] ==
                        Approx(shellset_ref[d][f12]).margin(abs_precision));
              }
            }
          };

          auto permute_deriv = [&](const int d) {
            assert(ncenters == 2);
            int dd;
            const auto xyz = d % 3;
            const auto atom = d / 3;

            // skip derivs w.r.t. dummy center
            assert(xsxs);
            if (atom == 1 || atom == 3) return -1;

            int a[2] = {0, 0};
            a[atom == 2 ? 1 : 0] = 1;
            int aa[2];
            assert(xsxs);
            std::tie(aa[0], aa[1]) = permuter(a[0], a[1]);
            if (aa[0])
              dd = xyz;
            else if (aa[1])
              dd = 3 + xyz;
            else
              assert(false);
            return dd;
          };

          // permuting geometric deriv index is order-specific
          int dd;
          switch (deriv_order) {
            case 0:
              dd = 0;
              break;

            case 1: {
              dd = permute_deriv(d);
              if (dd == -1) return;  // skip unneeded derivatives
            } break;

            case 2: {
              int d1, d2;
              std::tie(d1, d2) = split_deriv2(ncenters, d);
              auto dd1 = permute_deriv(d1);
              auto dd2 = permute_deriv(d2);
              if (dd1 == -1 || dd2 == -1) return;  // skip unneeded derivatives
              dd = merge_deriv2(ncenters, dd1, dd2);
            } break;

            default:
              assert(false);  // can't compare yet deriv_order > 2
          }
          compare(dd);

        }  // d
      };

      validate(permute2_12);
    }
  }
#endif  // LIBINT2_SUPPORT_ERI2
}

TEST_CASE_METHOD(libint2::unit::DefaultFixture,
                 "2-e 4-c integrals permute correctly", "[engine][permute]") {
  SECTION("deriv_order=0") { validate4<0>(obs, atoms); }  // section
  SECTION("deriv_order=1") { validate4<1>(obs, atoms); }  // section
  SECTION("deriv_order=2") { validate4<2>(obs, atoms); }  // section
}

TEST_CASE_METHOD(libint2::unit::DefaultFixture,
                 "2-e 3-c integrals permute correctly", "[engine][permute]") {
  SECTION("deriv_order=0") { validate3<0>(obs, dfbs, atoms); }  // section
  SECTION("deriv_order=1") { validate3<1>(obs, dfbs, atoms); }  // section
  SECTION("deriv_order=2") { validate3<2>(obs, dfbs, atoms); }  // section
}

TEST_CASE_METHOD(libint2::unit::DefaultFixture,
                 "2-e 2-c integrals permute correctly", "[engine][permute]") {
  SECTION("deriv_order=0") { validate2<0>(dfbs, atoms); }  // section
  SECTION("deriv_order=1") { validate2<1>(dfbs, atoms); }  // section
  SECTION("deriv_order=2") { validate2<2>(dfbs, atoms); }  // section
}
