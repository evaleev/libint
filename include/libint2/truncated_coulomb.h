/* Binary64 truncated-Coulomb radial evaluator.
 * SPDX-License-Identifier: LGPL-3.0-or-later
 *
 * No extended precision, MPFR, or high-precision fallback.
 */
#ifndef LIBINT2_TRUNCATED_COULOMB_H
#define LIBINT2_TRUNCATED_COULOMB_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace libint2 {
namespace detail {
namespace truncated_coulomb {

struct Workspace;  // defined in truncated_coulomb.cc

/// Mutable work storage of os_core_ints::truncated_coulomb_gm_eval: series
/// coefficients and recurrence buffers reused between calls. Not thread-safe;
/// Engine keeps one per Engine object.
class Scratch {
 public:
  explicit Scratch(int mmax);
  ~Scratch();
  Scratch(const Scratch&);
  Scratch& operator=(const Scratch&);
  Scratch(Scratch&&) noexcept;
  Scratch& operator=(Scratch&&) noexcept;

 private:
  std::unique_ptr<Workspace> workspace_;
  friend void evaluate(double*, double, double, int, double, Scratch&);
};

/// Computes \f$ G_m(R,T) \f$, m = 0..mmax, for the kernel
/// \f$ \theta(R - r)/r \f$ in units of \f$ \rho^{-1/2} \f$:
/// R = 0 gives zeros, R = infinity the Boys function \f$ F_m(T) \f$.
/// @param[out] Gm array of at least mmax+1 values
/// @param[in] precision target absolute error, >= double epsilon; the
///            F_m and zero limits are used only when exact to precision
///            relative to F_m
/// @throw std::invalid_argument for R < 0, NaN, or T < 0 or not finite
void evaluate(double* Gm, double R, double T, int mmax, double precision,
              Scratch& scratch);
}  // namespace truncated_coulomb
}  // namespace detail

namespace os_core_ints {
/// Interpolation-free core integral evaluator for the truncated Coulomb
/// kernel \f$ \theta(R_c - r_{12})/r_{12} \f$. All arithmetic is binary64.
/// The evaluator is immutable and can be shared between threads; mutable
/// storage is in a scratch_type object, one per thread.
template <typename Real>
class truncated_coulomb_gm_eval {
 public:
  typedef Real value_type;
  typedef detail::truncated_coulomb::Scratch scratch_type;

  /// @param precision target absolute precision; values below double
  ///        epsilon are raised to it
  truncated_coulomb_gm_eval(int mmax, Real precision)
      : mmax_(mmax),
        precision_(std::max(static_cast<double>(precision),
                            std::numeric_limits<double>::epsilon())) {
    if (mmax < 0 || !std::isfinite(precision_))
      throw std::invalid_argument("invalid truncated Coulomb evaluator setup");
  }

  static std::shared_ptr<const truncated_coulomb_gm_eval> instance(
      int mmax, Real precision = std::numeric_limits<Real>::epsilon()) {
    return std::make_shared<const truncated_coulomb_gm_eval>(mmax, precision);
  }

  /// @return the maximum m for which \f$ G_m \f$ can be computed
  int max_m() const { return mmax_; }
  /// @return the target absolute precision
  Real precision() const { return static_cast<Real>(precision_); }

  /// @param[out] Gm array of at least mmax+1 values
  /// @param[in] rho \f$ \rho = \zeta \eta / (\zeta + \eta) \f$
  /// @param[in] T \f$ \rho |PQ|^2 \f$
  /// @param[in] cutoff \f$ R_c \geq 0 \f$, in length units; infinity gives
  ///            Coulomb
  void eval(Real* Gm, Real rho, Real T, int mmax, Real cutoff,
            scratch_type& scratch) const {
    if (mmax < 0 || mmax > mmax_ || !(rho > 0) || !std::isfinite(rho) ||
        !(cutoff >= 0))
      throw std::invalid_argument(
          "invalid truncated Coulomb evaluator request");
    const double R =
        static_cast<double>(cutoff) * std::sqrt(static_cast<double>(rho));
    evaluate(Gm, R, static_cast<double>(T), mmax, scratch);
  }
  /// same as above, with temporary scratch (allocates on every call)
  void eval(Real* Gm, Real rho, Real T, int mmax, Real cutoff) const {
    scratch_type scratch(mmax);
    eval(Gm, rho, T, mmax, cutoff, scratch);
  }

 private:
  int mmax_;
  double precision_;

  void evaluate(double* Gm, double R, double T, int mmax,
                scratch_type& scratch) const {
    detail::truncated_coulomb::evaluate(Gm, R, T, mmax, precision_, scratch);
  }
  template <typename Output>
  void evaluate(Output* Gm, double R, double T, int mmax,
                scratch_type& scratch) const {
    std::vector<double> values(mmax + 1);
    detail::truncated_coulomb::evaluate(values.data(), R, T, mmax, precision_,
                                        scratch);
    std::copy(values.begin(), values.end(), Gm);
  }
};
}  // namespace os_core_ints
}  // namespace libint2
#endif
