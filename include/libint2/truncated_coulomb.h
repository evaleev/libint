/* Binary64 truncated-Coulomb radial evaluator.
 * SPDX-License-Identifier: LGPL-3.0-or-later
 *
 * No extended precision, MPFR, or high-precision fallback.
 */
#ifndef LIBINT2_TRUNCATED_COULOMB_H
#define LIBINT2_TRUNCATED_COULOMB_H

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

namespace libint2 {
namespace detail {
namespace truncated_coulomb {

struct Options {
  double precision = std::numeric_limits<double>::epsilon();
  unsigned max_extra_orders = 4096;
  unsigned max_taylor_terms = 4096;
};

struct Workspace;
struct Scratch {
  explicit Scratch(int mmax);
  ~Scratch();
  Scratch(const Scratch&);
  Scratch& operator=(const Scratch&);
  Scratch(Scratch&&) noexcept;
  Scratch& operator=(Scratch&&) noexcept;

 private:
  std::unique_ptr<Workspace> workspace_;
  friend const double* evaluate_into(double, double, int, Scratch&,
                                     const Options&);
  friend void evaluate_to(double*, double, double, int, Scratch&,
                          const Options&);
  friend void evaluate_to_unchecked(double*, double, double, int, Scratch&,
                                    const Options&);
};

// Production entry point. Storage grows only when the caller's workspace
// capacity is exceeded; the returned array lives until its next evaluation.
const double* evaluate_into(double R, double T, int mmax, Scratch& scratch,
                            const Options& options = Options());
// Write directly to caller-owned storage of at least mmax+1 doubles.
void evaluate_to(double* output, double R, double T, int mmax, Scratch& scratch,
                 const Options& options = Options());
// Engine-only entry: setup and shell data establish the checked API's
// preconditions; output has at least mmax+1 elements and scratch is nonnull.
void evaluate_to_unchecked(double* output, double R, double T, int mmax,
                           Scratch& scratch, const Options& options);
}  // namespace truncated_coulomb
}  // namespace detail

namespace os_core_ints {
/// Interpolation-free radial G_m for theta(Rc-r12)/r12.
/// All runtime numerical arithmetic uses double. Geometry and normalization
/// remain with the caller; the evaluator itself has no mutable work arrays.
template <typename Real>
struct truncated_coulomb_gm_eval {
  using value_type = Real;
  using scratch_type = detail::truncated_coulomb::Scratch;

  truncated_coulomb_gm_eval(int mmax, Real precision) : mmax_(mmax) {
    options_.precision = static_cast<double>(precision);
    if (mmax < 0 ||
        !(options_.precision >= std::numeric_limits<double>::epsilon()) ||
        !std::isfinite(options_.precision))
      throw std::invalid_argument("invalid truncated Coulomb evaluator setup");
  }
  void operator()(Real* output, Real rho, Real T, int mmax, Real cutoff,
                  scratch_type* scratch = nullptr) const {
    if (!output || mmax < 0 || mmax > mmax_ || !(rho > 0) ||
        !std::isfinite(rho) || !(T >= 0) || !std::isfinite(T) || !(cutoff >= 0))
      throw std::invalid_argument(
          "invalid truncated Coulomb evaluator request");
    const double R =
        static_cast<double>(cutoff) * std::sqrt(static_cast<double>(rho));
    if (!scratch) {
      scratch_type local(mmax);
      evaluate_output(output, R, static_cast<double>(T), mmax, local, options_);
      return;
    }
    evaluate_output(output, R, static_cast<double>(T), mmax, *scratch,
                    options_);
  }
  // Internal Engine path. Constructor and Engine parameter validation have
  // established the request contract, so no per-primitive scans are needed.
  void eval_unchecked(Real* output, Real rho, Real T, int mmax, Real cutoff,
                      scratch_type* scratch) const {
    evaluate_unchecked(
        output,
        static_cast<double>(cutoff) * std::sqrt(static_cast<double>(rho)),
        static_cast<double>(T), mmax, *scratch, options_);
  }

 private:
  static void evaluate_unchecked(
      double* output, double R, double T, int mmax, scratch_type& scratch,
      const detail::truncated_coulomb::Options& options) {
    detail::truncated_coulomb::evaluate_to_unchecked(output, R, T, mmax,
                                                     scratch, options);
  }
  template <typename Output>
  static void evaluate_unchecked(
      Output* output, double R, double T, int mmax, scratch_type& scratch,
      const detail::truncated_coulomb::Options& options) {
    const double* values =
        detail::truncated_coulomb::evaluate_into(R, T, mmax, scratch, options);
    for (int m = 0; m <= mmax; ++m) output[m] = static_cast<Output>(values[m]);
  }
  static void evaluate_output(
      double* output, double R, double T, int mmax, scratch_type& scratch,
      const detail::truncated_coulomb::Options& options) {
    detail::truncated_coulomb::evaluate_to(output, R, T, mmax, scratch,
                                           options);
    for (int m = 0; m <= mmax; ++m)
      if (!std::isfinite(output[m]))
        throw std::runtime_error("nonfinite truncated Coulomb result");
  }
  template <typename Output>
  static void evaluate_output(
      Output* output, double R, double T, int mmax, scratch_type& scratch,
      const detail::truncated_coulomb::Options& options) {
    const double* values =
        detail::truncated_coulomb::evaluate_into(R, T, mmax, scratch, options);
    for (int m = 0; m <= mmax; ++m) {
      if (!std::isfinite(values[m]))
        throw std::runtime_error("nonfinite truncated Coulomb result");
      output[m] = static_cast<Output>(values[m]);
    }
  }

  int mmax_;
  detail::truncated_coulomb::Options options_;
};
}  // namespace os_core_ints
}  // namespace libint2
#endif
