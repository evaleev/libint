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

#include <libint2/solidharmonics.h>

#include <array>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "catch.hpp"
#include "fixture.h"

TEST_CASE("SolidHarmonicsCoefficients ctor", "[shell]") {
  using libint2::solidharmonics::SolidHarmonicsCoefficients;
  REQUIRE_NOTHROW(SolidHarmonicsCoefficients<long double>::instance(10));
  auto sh12 = SolidHarmonicsCoefficients<long double>::instance(12);
  CHECK_NOTHROW(sh12.coeff(12, 0, 0, 0, 12));
  CHECK(sh12.coeff(12, 0, 0, 0, 12) == Approx(1));
  CHECK_NOTHROW(
      SolidHarmonicsCoefficients<long double>::make_instance(25).coeff(25, 0, 0,
                                                                       0, 25));
  CHECK(SolidHarmonicsCoefficients<long double>::make_instance(25).coeff(
            25, 0, 0, 0, 25) == Approx(1));
}

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
  REQUIRE(!s1.contr[0].pure);
  REQUIRE(s1.contr[0].coeff.size() == 1);
  REQUIRE(s1.contr[0].coeff[0] ==
          Approx(1.64592278064949));  // (2./\[Pi])^(3/4) 2^l
                                      // \[Alpha]^((2l+3)/4) / Sqrt[(2l-1)!!]
  REQUIRE(s1.O == std::array<double, 3>{0, 0, 0});
}

TEST_CASE("Engine ctor", "[engine]") {
  REQUIRE_NOTHROW(Engine{});
  auto a = Engine{};
}

namespace libint2 {
namespace detail {
struct EngineLayoutInspector {
  /// @return the offset of the first byte of @p member within @p engine
  template <typename T>
  static std::ptrdiff_t begin(const Engine& engine, const T& member) {
    return reinterpret_cast<const char*>(&member) -
           reinterpret_cast<const char*>(&engine);
  }
  /// @return the offset of the byte past @p member within @p engine
  template <typename T>
  static std::ptrdiff_t end(const Engine& engine, const T& member) {
    return begin(engine, member) + static_cast<std::ptrdiff_t>(sizeof(T));
  }

  // the state of the Engine, in particular spbra_ and spket_ that compute()
  // writes, must keep a distance of at least LIBINT2_ENGINE_FALSE_SHARING_GUARD
  // bytes from either end of the object, see
  // https://github.com/evaleev/libint/issues/478
  static void check(const Engine& e) {
    constexpr std::ptrdiff_t guard = LIBINT2_ENGINE_FALSE_SHARING_GUARD;
    constexpr std::ptrdiff_t size = sizeof(Engine);
    CHECK(begin(e, e.oper_) >= guard);
    CHECK(begin(e, e.spbra_) >= guard);
    CHECK(end(e, e.spket_) <= size - guard);
    CHECK(end(e, e.buildfnptrs_) <= size - guard);
  }
};
}  // namespace detail
}  // namespace libint2

TEST_CASE("Engine keeps its state away from its neighbors", "[engine]") {
  libint2::detail::EngineLayoutInspector::check(Engine{});
}

TEST_CASE("Engine::set", "[engine]") {
  REQUIRE_THROWS_AS(Engine{}.set(Operator::overlap),
                    Engine::using_default_initialized);
  REQUIRE_THROWS_AS(Engine{}.set(BraKet::x_x),
                    Engine::using_default_initialized);
  REQUIRE_NOTHROW(Engine{}.set(CartesianShellNormalization::uniform));
  REQUIRE_NOTHROW(Engine{}.prescale_by(1.5));
  REQUIRE_NOTHROW(Engine(Operator::overlap, 1, 0)
                      .set(CartesianShellNormalization::uniform)
                      .set_precision(1e-20)
                      .set(Operator::overlap)
                      .prescale_by(1.3)
                      .set(BraKet::x_x));
}

namespace {

// the task name of each Operator, in the order of Operator; an independent copy
// of LIBINT2_PP_NBODY_OPERATOR_LIST in engine.impl.h, from which the dispatch
// tables of Engine are generated
constexpr const char* task_names[] = {
    "overlap",       "kinetic",       "elecpot",
    "elecpot",       "elecpot",       "elecpot",
    "elecpot",       "1emultipole",   "2emultipole",
    "3emultipole",   "sphemultipole", "opVop",
    "opVop",         "opemuop",       "giao_overlap",
    "giao_kinetic",  "giao_elecpot",  "angmom",
    "eri",           "eri",           "eri",
    "eri",           "eri",           "eri",
    "eri",           "eri",           "eri",
    "eri",           "eri",           "eri",
    "coulomb_opop",  "coulomb_op",    "opop_coulomb_opop",
    "op_coulomb_op", "eri",           "giao_eri"};
static_assert(sizeof(task_names) / sizeof(task_names[0]) ==
                  static_cast<std::size_t>(Operator::last_oper) + 1,
              "add the task name of the new Operator to task_names");

/// @return the name of the task that computes the integrals of operator @p o
/// over braket @p bk with derivative order @p d, as reported by Engine (the
/// number of centers is omitted if it is the default for the rank of @p o)
std::string task_name(Operator o, BraKet bk, int d) {
  const auto default_ncenter = libint2::rank(o) == 1 ? 2 : 4;
  const auto ncenter = libint2::rank(bk);
  return (ncenter != default_ncenter ? std::to_string(ncenter)
                                     : std::string{}) +
         task_names[static_cast<int>(o)] +
         (d > 0 ? std::to_string(d) : std::string{});
}

// whether Engine::compute2_ptrs()[I] is Engine::compute2<o,bk,d> of the
// (o,bk,d) that I encodes, for each I
template <std::size_t... I>
std::array<bool, sizeof...(I)> compute2_ptrs_match(
    const std::vector<Engine::compute2_ptr_type>& ptrs,
    std::index_sequence<I...>) {
  constexpr auto o0 = static_cast<std::size_t>(Operator::first_2body_oper);
  constexpr auto bk0 = static_cast<std::size_t>(BraKet::first_2body_braket);
  constexpr auto nbk = libint2::nbrakets_2body;
  constexpr auto nd = libint2::nderivorders_2body;
  return {{(ptrs.at(I) ==
            static_cast<Engine::compute2_ptr_type>(
                &Engine::compute2<static_cast<Operator>(o0 + I / (nbk * nd)),
                                  static_cast<BraKet>(bk0 + (I / nd) % nbk),
                                  I % nd>))...}};
}

}  // namespace

TEST_CASE("Engine::compute2_ptrs", "[engine]") {
  constexpr auto n = libint2::nopers_2body * libint2::nbrakets_2body *
                     libint2::nderivorders_2body;
  const auto& ptrs = Engine{}.compute2_ptrs();
  REQUIRE(ptrs.size() == n);
  const auto match = compute2_ptrs_match(ptrs, std::make_index_sequence<n>{});
  for (std::size_t i = 0; i != n; ++i) {
    INFO("compute2_ptrs()[" << i << "]");
    CHECK(match[i]);
  }
}

// _initialize() selects the kernels for (Operator, rank of BraKet, derivative
// order). With max_l beyond the kernels it throws lmax_exceeded, which names
// the selected task. If the library has no kernels for the combination it
// throws nothing, or Engine rejects the combination (std::invalid_argument).
TEST_CASE("Engine selects the kernels of each operator", "[engine]") {
  // beyond the kernels of any task: LIBINT2_MAX_AM is not, e.g. the 3-center
  // ERIs can have higher angular momenta (LIBINT2_ERI3_MAX_AM); without kernels
  // Engine builds the core evaluators for max_l, and those with tables throw
  // std::invalid_argument
  const auto max_l = 64;
  std::set<std::string> tasks_found;
  for (int oi = static_cast<int>(Operator::first_oper);
       oi <= static_cast<int>(Operator::last_oper); ++oi) {
    const auto o = static_cast<Operator>(oi);
    const auto brakets =
        libint2::rank(o) == 1
            ? std::vector<BraKet>{BraKet::x_x}
            : std::vector<BraKet>{BraKet::xx_xx, BraKet::xs_xx, BraKet::xx_xs,
                                  BraKet::xs_xs};
    for (auto bk : brakets) {
      for (int d = 0; d <= LIBINT2_MAX_DERIV_ORDER; ++d) {
        const auto expected = task_name(o, bk, d);
        INFO("Operator #" << oi << " BraKet #" << static_cast<int>(bk)
                          << " deriv_order=" << d << ": expected task "
                          << expected);
        try {
          Engine(o, 1, max_l, d, std::numeric_limits<double>::epsilon(),
                 libint2::default_params(o), bk);
        } catch (const Engine::lmax_exceeded& e) {
          CHECK(std::string(e.task_name()) == expected);
          tasks_found.insert(e.task_name());
        } catch (const std::exception&) {
          // the combination is rejected, or there are no kernels for it
        }
      }
    }
  }
  // the kernels of these tasks, if the library has them, must be found
  // (LIBINT2_TASK_EXISTS_* always include the number of centers)
  const auto check_found = [&](const char* task) {
    INFO("task " << task);
    CHECK(tasks_found.count(task) == 1);
  };
#if LIBINT2_TASK_EXISTS_2overlap
  check_found("overlap");
#endif
#if LIBINT2_TASK_EXISTS_2elecpot1
  check_found("elecpot1");
#endif
#if LIBINT2_TASK_EXISTS_4eri
  check_found("eri");
#endif
#if LIBINT2_TASK_EXISTS_4eri1
  check_found("eri1");
#endif
#if LIBINT2_TASK_EXISTS_3eri
  check_found("3eri");
#endif
#if LIBINT2_TASK_EXISTS_2eri
  check_found("2eri");
#endif
}

template <typename RealBuf>
void check_uniform(int l, RealBuf&& S) {
  const auto n = (l + 1) * (l + 2) / 2;
  for (int i = 0; i != n; ++i) {
    REQUIRE(S[i * n + i] == Approx(1.));
  }
};

template <typename RealBuf>
void check_std_2(RealBuf&& S) {
  REQUIRE(S[0] == Approx(1.));
  REQUIRE(S[7] == Approx(1. / 3));
  REQUIRE(S[14] == Approx(1. / 3));
  REQUIRE(S[21] == Approx(1.));
  REQUIRE(S[28] == Approx(1. / 3));
  REQUIRE(S[35] == Approx(1.));
};

template <typename RealBuf>
void check_std_3(RealBuf&& S) {
  REQUIRE(S[0] == Approx(1.));
  REQUIRE(S[11] == Approx(1. / 5));
  REQUIRE(S[22] == Approx(1. / 5));
  REQUIRE(S[33] == Approx(1. / 5));
  REQUIRE(S[44] == Approx(1. / 15));
  REQUIRE(S[55] == Approx(1. / 5));
  REQUIRE(S[66] == Approx(1.));
  REQUIRE(S[77] == Approx(1. / 5));
  REQUIRE(S[88] == Approx(1. / 5));
  REQUIRE(S[99] == Approx(1.));
};

TEST_CASE("cartesian uniform normalization", "[engine][conventions]") {
#if defined(LIBINT2_SUPPORT_ONEBODY) && defined(LIBINT2_SUPPORT_ERI)
  if (LIBINT_CGSHELL_ORDERING != LIBINT_CGSHELL_ORDERING_STANDARD) return;

  std::vector<Shell> obs{Shell{{1.0}, {{2, false, {1.0}}}, {{0.0, 0.0, 0.0}}},
                         Shell{{1.0}, {{3, false, {1.0}}}, {{0.0, 0.0, 0.0}}}};
  {
    const auto lmax = std::min(3, LIBINT2_MAX_AM_overlap);
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
    const auto lmax = std::min(3, LIBINT2_MAX_AM_eri);
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
