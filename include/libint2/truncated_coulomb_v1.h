/* Independent truncated-Coulomb evaluator. LGPL-3.0-or-later. */
#ifndef LIBINT2_TRUNCATED_COULOMB_V1_H
#define LIBINT2_TRUNCATED_COULOMB_V1_H
#include <libint2/truncated_coulomb_table.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace libint2 {
namespace os_core_ints {
/// Radial G_m for theta(Rc-r)/r. Immutable, independently generated
/// Chebyshev patches cover the transition; the exterior uses limiting values.
template <typename Real>
struct truncated_coulomb_v1_gm_eval {
  typedef Real value_type;
  truncated_coulomb_v1_gm_eval(unsigned, Real precision) {
    if (!(precision >= std::numeric_limits<double>::epsilon()))
      throw std::invalid_argument(
          "truncated Coulomb precision below double epsilon");
  }
  void operator()(Real* result, Real rho, Real T, int mmax, Real cutoff) const {
    using std::isfinite;
    using std::sqrt;
    if (mmax < 0 || !(rho > 0) || !isfinite(rho) || !(T >= 0) || !isfinite(T) ||
        !(cutoff >= 0))
      throw std::invalid_argument("invalid truncated Coulomb arguments");
    if (cutoff == 0) {
      std::fill(result, result + mmax + 1, Real(0));
      return;
    }
    const long double t = static_cast<long double>(T);
    if (!isfinite(cutoff)) {
      boys(result, t, mmax);
      return;
    }
    if (mmax > detail::truncated_coulomb_max_m)
      throw std::invalid_argument(
          "truncated Coulomb order exceeds generated table");
    const long double r =
        static_cast<long double>(cutoff) * sqrt(static_cast<long double>(rho));
    const long double s = sqrt(t), w = detail::truncated_coulomb_width;
    if (r - s > w) {
      boys(result, t, mmax);
      return;
    }
    if (s - r > w) {
      std::fill(result, result + mmax + 1, Real(0));
      return;
    }
    const bool far = r > 8;
    const long double u = far ? 8 / r : r / 8;
    const long double v = far ? (s - r + w) / (2 * w) : s / (8 + w);
    int index = far ? 1 : 0;
    while (detail::truncated_coulomb_nodes[index].axis >= 0) {
      const auto& node = detail::truncated_coulomb_nodes[index];
      index = (node.axis == 0 ? u : v) <= node.split ? node.left : node.right;
    }
    const auto& leaf = detail::truncated_coulomb_nodes[index];
    const long double x = (2 * u - leaf.a - leaf.b) / (leaf.b - leaf.a);
    const long double y = (2 * v - leaf.c - leaf.d) / (leaf.d - leaf.c);
    const long double factor = far ? u / 8 : 1;
    for (int m = 0; m <= mmax; ++m) {
      const auto& polynomial =
          detail::truncated_coulomb_polynomials[leaf.poly + m];
      if (polynomial.nx == 0) {
        result[m] = Real(0);
        continue;
      }
      const double* c =
          detail::truncated_coulomb_coefficients
              [polynomial.offset / detail::truncated_coulomb_block_size] +
          polynomial.offset % detail::truncated_coulomb_block_size;
      long double rows[detail::truncated_coulomb_degree + 1];
      for (unsigned i = 0; i < polynomial.nx; ++i) {
        long double b1 = 0, b2 = 0;
        for (int j = polynomial.ny[i] - 1; j > 0; --j) {
          const auto b = 2 * y * b1 - b2 + c[j];
          b2 = b1;
          b1 = b;
        }
        rows[i] = polynomial.ny[i] ? y * b1 - b2 + c[0] : 0;
        c += polynomial.ny[i];
      }
      long double b1 = 0, b2 = 0;
      for (int i = polynomial.nx - 1; i > 0; --i) {
        const auto b = 2 * x * b1 - b2 + rows[i];
        b2 = b1;
        b1 = b;
      }
      result[m] = static_cast<Real>((x * b1 - b2 + rows[0]) * factor);
    }
  }

 private:
  static void boys(Real* result, long double t, int mmax) {
    using std::erf;
    using std::exp;
    using std::sqrt;
    if (t <= mmax / 2.L + 4) {
      // Positive-term incomplete-gamma series at the highest order, followed
      // by stable downward recurrence. This has no tabulated m limit.
      long double term = 1 / (2.L * mmax + 1), sum = term;
      for (int k = 1;; ++k) {
        term *= t / (mmax + k + 0.5L);
        sum += term;
        if (term <= sum * std::numeric_limits<long double>::epsilon()) break;
      }
      long double g = exp(-t) * sum;
      const auto et = exp(-t);
      result[mmax] = static_cast<Real>(g);
      for (int m = mmax - 1; m >= 0; --m) {
        g = (2 * t * g + et) / (2 * m + 1);
        result[m] = static_cast<Real>(g);
      }
    } else {
      const long double pi = acos(-1.L);
      long double g = sqrt(pi) * erf(sqrt(t)) / (2 * sqrt(t));
      const auto et = exp(-t);
      result[0] = static_cast<Real>(g);
      for (int m = 0; m < mmax; ++m) {
        g = ((2 * m + 1) * g - et) / (2 * t);
        result[m + 1] = static_cast<Real>(g);
      }
    }
  }
};
}  // namespace os_core_ints
}  // namespace libint2
#endif
