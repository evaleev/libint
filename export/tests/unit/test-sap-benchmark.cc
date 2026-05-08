// SAP 2c (q_gau) vs 3c (Coulomb xs_xx · SAP-S) benchmark.
// Tagged [.benchmark] so the default unit_tests-libint2 run skips it.
// Opt-in: ./unit_tests-libint2 "[.benchmark]"
//
// Boys / RR / tform / q_gau-total times are read directly from the
// thread-local timers wired by LIBINT2_BOYS_TIMING (see
// include/libint2/util/timer.h). The 1-body path uses the
// `eng.compute1(s1, s2, &sp)` overload added on this branch so the per-shell
// pair geometry that production code amortises via ShellPair is amortised
// here too — the "rest" bucket reflects production usage rather than
// per-call recomputation.

#include "catch.hpp"

#include <libint2.hpp>
#include <libint2/util/timer.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

using libint2::Atom;
using libint2::BraKet;
using libint2::Engine;
using libint2::GaussianPotentialCentersData;
using libint2::GaussianPotentialData;
using libint2::NuclearModel;
using libint2::Operator;
using libint2::ScreeningMethod;
using libint2::Shell;
using libint2::ShellPair;

constexpr int Lmax_obs = 5;  // S..H

// Fe AHGBSP3-9 from BSE. All shells are uncontracted single primitives.
// 35 S, 28 P, 26 D, 26 F, 26 G, 26 H = 167 total.
const std::array<std::vector<double>, Lmax_obs + 1>& fe_ahgbsp3_9_alphas() {
  static const std::array<std::vector<double>, Lmax_obs + 1> data = {{
      // S (35)
      {5.1525236671e+07, 2.6714040435e+07, 1.3850299435e+07, 7.1808978097e+06,
       3.7230453821e+06, 1.9302693458e+06, 1.0007774187e+06, 5.1886823151e+05,
       2.6901510430e+05, 1.3947496098e+05, 7.2312908934e+04, 3.7491724406e+04,
       1.9438153155e+04, 1.0078005322e+04, 5.2250947122e+03, 2.7090296025e+03,
       1.4045374853e+03, 7.2820376185e+02, 3.7754828500e+02, 1.9574563463e+02,
       1.0148729315e+02, 5.2617626390e+01, 2.7280406453e+01, 1.4143940487e+01,
       7.3331404668e+00, 3.8019778967e+00, 1.9711931051e+00, 1.0219949624e+00,
       5.2986878887e-01, 2.7471851011e-01, 1.4243197822e-01, 7.3846019374e-02,
       3.8286588767e-02, 1.9850262639e-02, 1.0291669734e-02},
      // P (28)
      {3.3644076128e+04, 1.8197365526e+04, 9.8425681487e+03, 5.3236358649e+03,
       2.8794414623e+03, 1.5574286719e+03, 8.4238005868e+02, 4.5562546528e+02,
       2.4643812786e+02, 1.3329314425e+02, 7.2095427996e+01, 3.8994884300e+01,
       2.1091503912e+01, 1.1407946074e+01, 6.1703155059e+00, 3.3373924802e+00,
       1.8051246417e+00, 9.7635354290e-01, 5.2808887471e-01, 2.8563204550e-01,
       1.5449230106e-01, 8.3561601234e-02, 4.5196693642e-02, 2.4445930738e-02,
       1.3222284231e-02, 7.1516524435e-03, 3.8681767671e-03, 2.0922145784e-03},
      // D (26)
      {1.3432353325e+03, 7.5727686018e+02, 4.2693058251e+02, 2.4069099674e+02,
       1.3569455617e+02, 7.6500628710e+01, 4.3128820772e+01, 2.4314769860e+01,
       1.3707957295e+01, 7.7281460728e+00, 4.3569031065e+00, 2.4562947569e+00,
       1.3847872641e+00, 7.8070262597e-01, 4.4013734527e-01, 2.4813658396e-01,
       1.3989216085e-01, 7.8867115664e-02, 4.4462977019e-02, 2.5066928197e-02,
       1.4132002204e-02, 7.9672102108e-03, 4.4916804872e-03, 2.5322783089e-03,
       1.4276245721e-03, 8.0485304944e-04},
      // F (26)
      {2.2630197192e+02, 1.3394148482e+02, 7.9276027533e+01, 4.6921150305e+01,
       2.7771249575e+01, 1.6436986262e+01, 9.7285689883e+00, 5.7580540039e+00,
       3.4080229016e+00, 2.0171085735e+00, 1.1938672699e+00, 7.0661494223e-01,
       4.1822461271e-01, 2.4753485417e-01, 1.4650860367e-01, 8.6714135755e-02,
       5.1323547912e-02, 3.0376899306e-02, 1.7979193742e-02, 1.0641356261e-02,
       6.2983059585e-03, 3.7277821524e-03, 2.2063646744e-03, 1.3058823927e-03,
       7.7291340060e-04, 4.5746472131e-04},
      // G (26)
      {6.0353965486e+01, 3.6744996859e+01, 2.2371268951e+01, 1.3620185529e+01,
       8.2923080607e+00, 5.0485636063e+00, 3.0736912208e+00, 1.8713397429e+00,
       1.1393182274e+00, 6.9364530318e-01, 4.2230853070e-01, 2.5711194798e-01,
       1.5653615541e-01, 9.5303108793e-02, 5.8022905456e-02, 3.5325789475e-02,
       2.1507220161e-02, 1.3094131113e-02, 7.9720330345e-03, 4.8535721962e-03,
       2.9549755955e-03, 1.7990627145e-03, 1.0953141730e-03, 6.6685453925e-04,
       4.0599764660e-04, 2.4718147563e-04},
      // H (26)
      {2.1199803458e+01, 1.3140857680e+01, 8.1454595045e+00, 5.0490243601e+00,
       3.1296757383e+00, 1.9399530540e+00, 1.2024944967e+00, 7.4537526128e-01,
       4.6202646389e-01, 2.8639057991e-01, 1.7752135575e-01, 1.1003794803e-01,
       6.8207850008e-02, 4.2279149021e-02, 2.6207048628e-02, 1.6244636273e-02,
       1.0069360018e-02, 6.2415685687e-03, 3.8688832387e-03, 2.3981563849e-03,
       1.4865152788e-03, 9.2142768007e-04, 5.7115388029e-04, 3.5403403004e-04,
       2.1945065726e-04, 1.3602814104e-04},
  }};
  return data;
}

std::size_t max_nprim_qgau(const GaussianPotentialCentersData& data) {
  std::size_t n = 0;
  for (const auto& ptr : data)
    if (ptr) n = std::max(n, ptr->size());
  return n;
}

// Measure the per-call cost of std::chrono::high_resolution_clock::now() on
// the current machine, in nanoseconds. Mirrors src/bin/profile/chrono.cc.
// Run once at test start so each per-component sub-run can subtract the
// right per-stop overhead — independent of host hardware.
std::size_t measure_chrono_now_overhead_ns() {
  using clock_t = std::chrono::high_resolution_clock;
  constexpr std::size_t nrepeats = 2'000'000;
  const auto t0 = clock_t::now();
  // Side-effect on a volatile sink so the compiler can't elide now().
  volatile clock_t::rep sink = 0;
  for (std::size_t i = 0; i != nrepeats; ++i) {
    sink ^= clock_t::now().time_since_epoch().count();
  }
  const auto t1 = clock_t::now();
  (void)sink;
  const double ns =
      std::chrono::duration<double, std::nano>(t1 - t0).count() / nrepeats;
  return static_cast<std::size_t>(ns + 0.5);
}

Shell make_shell(double alpha, int L, std::array<double, 3> O) {
  libint2::svector<double> a{alpha};
  libint2::svector<double> c{1.0};
  return Shell{a, {{L, /*pure=*/true, c}}, O};
}

}  // namespace

TEST_CASE("SAP 2c-vs-3c benchmark on Fe AHGBSP3-9",
          "[.benchmark][engine][sap]") {
  // Re-measure chrono::now() overhead on this host so component timers
  // subtract the right per-stop cost (calibrates ~16 ns on Apple M2,
  // ~30–40 ns on x86 with vDSO clock). 2 M reps takes ~30 ms.
  const std::size_t chrono_overhead_ns = measure_chrono_now_overhead_ns();
  std::cerr << "Measured chrono::now() overhead = " << chrono_overhead_ns
            << " ns (used for per-stop overhead correction)\n";

  constexpr int Z_Fe = 26;
  const std::vector<Atom> atoms{{Z_Fe, 0.0, 0.0, 0.0}};
  const std::array<double, 3> nu_pos{{0.5, 0.3, 0.2}};

  // 2c side: SAP-only q_gau data (zero out the point-charge slot so the
  // operator is purely Gaussian potential).
  auto q_gau_full =
      libint2::make_q_gau_data(NuclearModel::PointCharge, atoms,
                               "sap_grasp_large");
  GaussianPotentialCentersData sap_only;
  sap_only.reserve(q_gau_full.size());
  for (const auto& ptr : q_gau_full) {
    if (ptr && !ptr->empty()) {
      auto data = *ptr;
      data[0].coefficient = 0.0;
      sap_only.push_back(
          std::make_shared<const GaussianPotentialData>(std::move(data)));
    } else {
      sap_only.push_back(ptr);
    }
  }
  auto point_charges = libint2::make_point_charges(atoms);

  // 3c side: SAP S-aux shell with prefactor adjustment so per-primitive
  // c_adj·(2α/π)^{3/4} = c_BSE·(α/π)^{3/2}.
  auto sap_by_Z = libint2::read_sap_basis_library(
      libint2::basis_data_path() + "/sap_grasp_large.g94");
  const auto& Fe_sap = sap_by_Z[Z_Fe];
  libint2::svector<double> aux_alphas, aux_coeffs;
  aux_alphas.reserve(Fe_sap.size());
  aux_coeffs.reserve(Fe_sap.size());
  const double inv_2_to_3_4 = std::pow(2.0, -0.75);
  for (const auto& p : Fe_sap) {
    aux_alphas.push_back(p.exponent);
    aux_coeffs.push_back(p.coefficient * std::pow(p.exponent / M_PI, 0.75) *
                         inv_2_to_3_4);
  }
  Shell::do_enforce_unit_normalization(false);
  Shell aux{aux_alphas, {{0, /*pure=*/true, aux_coeffs}}, {{0.0, 0.0, 0.0}}};
  Shell::do_enforce_unit_normalization(true);

  const auto& alphas_by_L = fe_ahgbsp3_9_alphas();

  // Engines.
  Engine eng2c(Operator::q_gau,
               std::max<std::size_t>(1, max_nprim_qgau(sap_only)), Lmax_obs);
  eng2c.set_params(std::make_tuple(sap_only, point_charges));

  Engine eng3c(Operator::coulomb,
               std::max<std::size_t>(aux.nprim(), 1), Lmax_obs, /*deriv=*/0);
  eng3c.set(BraKet::xs_xx);

  Engine eng3c_noscreen(Operator::coulomb,
                        std::max<std::size_t>(aux.nprim(), 1), Lmax_obs,
                        /*deriv=*/0);
  eng3c_noscreen.set(BraKet::xs_xx);
  eng3c_noscreen.set_precision(0.0);

  // Per-cell measurement.
  // side: 0=2c (q_gau), 1=3c default screening, 2=3c no screening.
  // Returns mean (Total, Boys, RR, tform) per shell-pair eval, each
  // measured in its OWN dedicated K-iteration sweep with only its timer
  // active — no nested chrono overhead between component timers. Total
  // is wall time of compute() with all inner timers disabled.
  auto bench_pair = [&](const std::vector<Shell>& mus,
                        const std::vector<Shell>& nus,
                        int side) -> std::array<double, 4> {
    Engine& eng = (side == 0) ? eng2c : (side == 1) ? eng3c : eng3c_noscreen;
    const bool is_2c = (side == 0);

    // Precompute ShellPair data once per (mu, nu). Same precomputed
    // ShellPair is fed into both compute1 (2c) and compute2 (3c).
    const double ln_prec = (eng.precision() <= 0.0)
                               ? std::numeric_limits<double>::lowest()
                               : std::log(eng.precision());
    std::vector<ShellPair> shellpairs;
    shellpairs.reserve(mus.size() * nus.size());
    for (const auto& mu : mus) {
      for (const auto& nu : nus) {
        shellpairs.emplace_back();
        shellpairs.back().init(mu, nu, ln_prec, ScreeningMethod::Original);
      }
    }

    auto do_one_compute = [&](std::size_t pp_idx, const Shell& mu,
                              const Shell& nu) {
      if (is_2c)
        eng.compute1(mu, nu, &shellpairs[pp_idx]);
      else
        eng.compute2<Operator::coulomb, BraKet::xs_xx, 0>(
            aux, Shell::unit(), mu, nu, /*spbra=*/nullptr,
            /*spket=*/&shellpairs[pp_idx]);
    };

    // Warmup pass with all timers off.
#ifdef LIBINT2_BOYS_TIMING
    ::libint2::detail::active_timer_kind() =
        ::libint2::detail::TimerKind::None;
#endif
    {
      double sink = 0.0;
      std::size_t pp = 0;
      for (const auto& mu : mus)
        for (const auto& nu : nus) {
          do_one_compute(pp, mu, nu);
          ++pp;
          const auto* buf = eng.results()[0];
          if (buf) sink += buf[0];
        }
      if (sink == 1.7976931348623157e+308) std::abort();
    }

    const std::size_t n_pairs = mus.size() * nus.size();
    const double target_s = 0.2;

    // Single sub-run: set the requested timer kind, do K passes, return
    // the measured value (in seconds). For Total (None) the wall clock is
    // returned; for Boys / RR / tform the matching timer's read.
#ifdef LIBINT2_BOYS_TIMING
    auto sub_run = [&](::libint2::detail::TimerKind kind, std::size_t K)
        -> double {
      double sink = 0.0;
      ::libint2::detail::boys_fm_timer().clear();
      ::libint2::detail::rr_kernel_timer().clear();
      ::libint2::detail::tform_timer().clear();
      // Use the runtime-measured chrono overhead so each stop() subtracts
      // the right per-call bias on this machine.
      ::libint2::detail::boys_fm_timer().set_now_overhead(chrono_overhead_ns);
      ::libint2::detail::rr_kernel_timer().set_now_overhead(chrono_overhead_ns);
      ::libint2::detail::tform_timer().set_now_overhead(chrono_overhead_ns);
      ::libint2::detail::active_timer_kind() = kind;
      const auto t0 = std::chrono::steady_clock::now();
      for (std::size_t k = 0; k < K; ++k) {
        std::size_t pp = 0;
        for (const auto& mu : mus)
          for (const auto& nu : nus) {
            do_one_compute(pp, mu, nu);
            ++pp;
            const auto* buf = eng.results()[0];
            if (buf) sink += buf[0];
          }
      }
      const auto t1 = std::chrono::steady_clock::now();
      ::libint2::detail::active_timer_kind() =
          ::libint2::detail::TimerKind::None;
      if (sink == 1.7976931348623157e+308) std::abort();
      switch (kind) {
        case ::libint2::detail::TimerKind::None:
          return std::chrono::duration<double>(t1 - t0).count();
        case ::libint2::detail::TimerKind::Boys:
          return ::libint2::detail::boys_fm_timer().read(0);
        case ::libint2::detail::TimerKind::RR:
          return ::libint2::detail::rr_kernel_timer().read(0);
        case ::libint2::detail::TimerKind::Tform:
          return ::libint2::detail::tform_timer().read(0);
      }
      return 0.0;
    };
#else
    auto sub_run = [&](int /*kind*/, std::size_t K) -> double {
      double sink = 0.0;
      const auto t0 = std::chrono::steady_clock::now();
      for (std::size_t k = 0; k < K; ++k) {
        std::size_t pp = 0;
        for (const auto& mu : mus)
          for (const auto& nu : nus) {
            do_one_compute(pp, mu, nu);
            ++pp;
            const auto* buf = eng.results()[0];
            if (buf) sink += buf[0];
          }
      }
      const auto t1 = std::chrono::steady_clock::now();
      if (sink == 1.7976931348623157e+308) std::abort();
      return std::chrono::duration<double>(t1 - t0).count();
    };
#endif

    // Probe with no inner timers to determine K from the total wall time.
#ifdef LIBINT2_BOYS_TIMING
    const auto k_total = ::libint2::detail::TimerKind::None;
    const auto k_boys = ::libint2::detail::TimerKind::Boys;
    const auto k_rr = ::libint2::detail::TimerKind::RR;
    const auto k_tform = ::libint2::detail::TimerKind::Tform;
#else
    const int k_total = 0, k_boys = 1, k_rr = 2, k_tform = 3;
#endif
    const double probe = sub_run(k_total, 1);
    const std::size_t K = std::max<std::size_t>(
        1, static_cast<std::size_t>(target_s / std::max(probe, 1e-9)));

    const double dt_total = (K > 1 ? sub_run(k_total, K) : probe);
    const double dt_boys = sub_run(k_boys, K);
    const double dt_rr = sub_run(k_rr, K);
    const double dt_tform = sub_run(k_tform, K);

    const double denom = static_cast<double>(K * n_pairs);
    return std::array<double, 4>{dt_total / denom, dt_boys / denom,
                                  dt_rr / denom, dt_tform / denom};
  };

  constexpr int N_REPS = 15;
  const char L_letter[] = {'s', 'p', 'd', 'f', 'g', 'h', 'i'};

  struct Cell {
    int L1, L2;
    double m2_t, m2_b, m2_r, m2_f;
    double m3_t, m3_b, m3_r, m3_f;
    double m3ns_t, m3ns_b, m3ns_r, m3ns_f;
    double speedup_total;  // 3c-scr   / 2c
    double speedup_noscr;  // 3c-noscr / 2c
  };
  std::vector<Cell> cells;

  for (int L1 = 0; L1 <= Lmax_obs; ++L1) {
    std::vector<Shell> mus;
    mus.reserve(alphas_by_L[L1].size());
    for (double a : alphas_by_L[L1])
      mus.push_back(make_shell(a, L1, {{0, 0, 0}}));

    for (int L2 = L1; L2 <= Lmax_obs; ++L2) {
      std::vector<Shell> nus;
      nus.reserve(alphas_by_L[L2].size());
      for (double a : alphas_by_L[L2])
        nus.push_back(make_shell(a, L2, nu_pos));

      std::array<std::array<std::vector<double>, 4>, 3> samples;
      for (auto& side : samples)
        for (auto& v : side) v.reserve(N_REPS);
      for (int rep = 0; rep < N_REPS; ++rep) {
        for (int s = 0; s < 3; ++s) {
          auto x = bench_pair(mus, nus, /*side=*/s);
          for (int i = 0; i < 4; ++i) samples[s][i].push_back(x[i]);
        }
      }
      auto mean = [](const std::vector<double>& v) {
        double s = 0;
        for (double x : v) s += x;
        return s / v.size();
      };
      double m[3][4];
      for (int s = 0; s < 3; ++s)
        for (int i = 0; i < 4; ++i) m[s][i] = mean(samples[s][i]) * 1e9;
      Cell c{L1,      L2,      m[0][0], m[0][1], m[0][2], m[0][3],
             m[1][0], m[1][1], m[1][2], m[1][3],
             m[2][0], m[2][1], m[2][2], m[2][3],
             m[1][0] / m[0][0], m[2][0] / m[0][0]};
      cells.push_back(c);
      std::cerr << "  (" << L_letter[L1] << "|V|" << L_letter[L2]
                << ") 2c: t=" << c.m2_t << " boys=" << c.m2_b
                << " rr=" << c.m2_r << " tform=" << c.m2_f
                << " | 3c: t=" << c.m3_t << " boys=" << c.m3_b
                << " rr=" << c.m3_r << " | 3c-noscr: t=" << c.m3ns_t
                << " | scr×=" << c.speedup_total
                << " noscr×=" << c.speedup_noscr << "\n";
    }
  }

  std::cout << "\nPer-integral timings (mean over " << N_REPS
            << " runs; ns per shell-pair eval).\n"
               "Three sides: 2c (q_gau), 3c (default screening), 3c-noscreen "
               "(set_precision(0.0)).\n"
               "Each component is measured in its OWN dedicated K-iteration "
               "sub-run with only that timer active (no nested chrono "
               "overhead). Total = wall time of compute() with all inner "
               "timers disabled. set_now_overhead = "
            << chrono_overhead_ns
            << " ns (measured at test startup on this host).\n";

  auto print_table = [&](const char* title, double Cell::*t,
                         double Cell::*b, double Cell::*r,
                         double Cell::*f) {
    std::cout << title << "\n";
    std::cout << std::string(64, '-') << "\n";
    std::cout << std::left << std::setw(18) << "Integral" << std::right
              << std::setw(10) << "Total" << std::setw(10) << "Boys"
              << std::setw(10) << "RR" << std::setw(10) << "tform" << "\n";
    std::cout << std::string(64, '-') << "\n";
    for (const auto& c : cells) {
      std::ostringstream label;
      label << "(" << L_letter[c.L1] << "|V_SAP|" << L_letter[c.L2] << ")";
      std::cout << std::left << std::setw(18) << label.str() << std::right
                << std::fixed << std::setprecision(1)
                << std::setw(10) << c.*t
                << std::setw(10) << c.*b
                << std::setw(10) << c.*r
                << std::setw(10) << c.*f << "\n";
    }
    std::cout << std::string(64, '-') << "\n\n";
  };
  print_table("=== 2c (q_gau) ===", &Cell::m2_t, &Cell::m2_b, &Cell::m2_r,
              &Cell::m2_f);
  print_table("=== 3c (coulomb xs_xx, default screening) ===", &Cell::m3_t,
              &Cell::m3_b, &Cell::m3_r, &Cell::m3_f);
  print_table("=== 3c (coulomb xs_xx, NO screening) ===", &Cell::m3ns_t,
              &Cell::m3ns_b, &Cell::m3ns_r, &Cell::m3ns_f);

  std::cout << "Speedups (3c_total / 2c_total) per cell:\n";
  std::cout << std::left << std::setw(20) << "Integral" << std::right
            << std::setw(14) << "scr (3c/2c)" << std::setw(14)
            << "no-scr (3c/2c)" << "\n";
  for (const auto& c : cells) {
    std::ostringstream label;
    label << "(" << L_letter[c.L1] << "|V_SAP|" << L_letter[c.L2] << ")";
    std::cout << std::left << std::setw(20) << label.str() << std::right
              << std::fixed << std::setprecision(2) << std::setw(13)
              << c.speedup_total << "x" << std::setw(13) << c.speedup_noscr
              << "x" << "\n";
  }
}
