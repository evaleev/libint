/* Binary64 truncated-Coulomb numerical kernel.
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#include <libint2/boys.h>
#include <libint2/truncated_coulomb.h>

// The error bounds below rely on IEEE semantics: infinities as sentinels,
// compensated sums, exact products and directed outward rounding. The build
// compiles this file with value-safe options (see CMakeLists.txt.export);
// fail here rather than return wrong integrals if they were overridden.
#if defined(__FAST_MATH__) || \
    (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__)
#error "truncated_coulomb.cc must be compiled without -ffast-math"
#endif

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <cstring>
#include <functional>
#include <vector>

namespace libint2 {
namespace detail {
namespace truncated_coulomb {
namespace {
constexpr double unit = std::numeric_limits<double>::epsilon() / 2;
// Target precision and resource limits of a request; exhausting a limit
// throws std::runtime_error.
struct Options {
  double precision;
  unsigned max_extra_orders = 4096;
  unsigned max_taylor_terms = 4096;
};
constexpr double infinity = std::numeric_limits<double>::infinity();
// Fused multiply-add is used only where the compiler emits it as an
// instruction. Elsewhere std::fma is a library call, emulated in software on
// CPUs without FMA (~150 ns on Ivy Bridge), so each use falls back to plain
// arithmetic with the extra rounding accounted for. Define
// LIBINT2_TRUNCATED_COULOMB_FAST_FMA=0 to test the fallback on FMA hardware.
#ifndef LIBINT2_TRUNCATED_COULOMB_FAST_FMA
#ifdef FP_FAST_FMA
#define LIBINT2_TRUNCATED_COULOMB_FAST_FMA 1
#else
#define LIBINT2_TRUNCATED_COULOMB_FAST_FMA 0
#endif
#endif
constexpr bool fast_fma = LIBINT2_TRUNCATED_COULOMB_FAST_FMA;
constexpr unsigned taylor_roundings = fast_fma ? 1 : 2;
// a*b+c, rounded once with fast FMA and twice otherwise; for arithmetic whose
// error budget does not depend on the number of roundings.
inline double mul_add(double a, double b, double c) {
  return fast_fma ? std::fma(a, b, c) : a * b + c;
}
// Preserve independent vector lanes on portable targets as well. A scalar
// libm fma call would serialize them; the budget accounts for both roundings.
inline double taylor_add(double weight, double origin, double sum) {
  return mul_add(weight, origin, sum);
}
// Numeric outward steps for nonnegative binary64 bounds. Exceptional and
// subnormal results retain nextafter's floating-point exception behavior.
inline double outward(double x, double direction) {
  static_assert(sizeof(double) == sizeof(std::uint64_t) &&
                    std::numeric_limits<double>::is_iec559 &&
                    std::numeric_limits<double>::digits == 53,
                "truncated Coulomb requires IEEE binary64 arithmetic");
  std::uint64_t bits;
  std::memcpy(&bits, &x, sizeof(bits));
  // Integer classification also leaves quiet NaNs quiet until nextafter.
  constexpr std::uint64_t min_normal = UINT64_C(0x0010000000000000);
  constexpr std::uint64_t max_finite = UINT64_C(0x7fefffffffffffff);
  if (bits >= min_normal && bits < max_finite &&
      (direction == infinity || bits > min_normal)) {
    if (direction == infinity)
      ++bits;
    else
      --bits;
    std::memcpy(&x, &bits, sizeof(x));
    return x;
  }
  return std::nextafter(x, direction);
}
// Bound on a*b+c for a*b, c >= 0, rounded toward direction (0 or infinity):
// one outward step per rounding.
inline double outward_mul_add(double a, double b, double c, double direction) {
  if (fast_fma) return outward(std::fma(a, b, c), direction);
  return outward(outward(a * b, direction) + c, direction);
}
// T - R^2 with the product formed exactly. Near R = sqrt(T) the result is the
// small displacement sqrt(T) - R times sqrt(T) + R, so rounding R^2 first
// would cost a relative error of order T / (T - R^2).
inline double difference_of_square(double T, double R) {
  if (fast_fma) return std::fma(-R, R, T);
  // Dekker's exact product R*R = p + e via Veltkamp splitting of R.
  const double p = R * R, c = 134217729. * R, hi = c - (c - R), lo = R - hi;
  const double e = ((hi * hi - p) + 2 * hi * lo) + lo * lo;
  return (T - p) - e;
}
// Double-word arithmetic: two_sum and two_product are error-free (Knuth,
// Dekker). Pairs are left unnormalized; callers bound their roundings.
struct Pair {
  double hi, lo;
};
inline Pair two_sum(double a, double b) {
  const double s = a + b, t = s - a;
  return {s, (a - (s - t)) + (b - t)};
}
inline Pair fast_two_sum(double a, double b) {  // |a| >= |b|
  const double s = a + b;
  return {s, b - (s - a)};
}
// Veltkamp splitting, exact for |a| < 2^996.
inline Pair split(double a) {
  const double c = 134217729. * a, hi = c - (c - a);
  return {hi, a - hi};
}
// a*b - p for p = fl(a*b), from split operands.
inline double product_error(Pair a, Pair b, double p) {
  return ((a.hi * b.hi - p) + a.hi * b.lo + a.lo * b.hi) + a.lo * b.lo;
}
inline Pair two_product(double a, double b) {
  const double p = a * b;
  if (fast_fma) return {p, std::fma(a, b, -p)};
  return {p, product_error(split(a), split(b), p)};
}
struct Sum {
  double value = 0, correction = 0, magnitude = 0;
  void add(double term) {
    const double adjusted = term - correction;
    const double next = value + adjusted;
    correction = (next - value) - adjusted;
    value = next;
    magnitude += std::abs(term);
  }
};
void validate(double R, double T, int mmax, const Options& options) {
  if (!(R >= 0) || !(T >= 0) || !std::isfinite(T) || mmax < 0 ||
      !(options.precision >= std::numeric_limits<double>::epsilon()) ||
      !std::isfinite(options.precision))
    throw std::invalid_argument("invalid truncated Coulomb request");
}
double root_difference(double R, double T, double s) {
  if (R >= .5 * s && R <= 2 * s) return difference_of_square(T, R) / (s + R);
  return s - R;
}
struct TaylorPlan {
  int degree;
  double remainder;
};
struct Candidate {
  double value;
  double error;
};
int zero_seed(double T, int n, double tolerance, unsigned limit) {
  double product = 1;
  const unsigned extra = std::min(limit, static_cast<unsigned>(INT_MAX - n));
  for (unsigned j = 1; j <= extra; ++j) {
    const int M = n + static_cast<int>(j);
    product = outward(product * outward(T / (M - .5), infinity), infinity);
    const double ratio = outward(T / (M + 1.5), infinity);
    const double denominator = outward(1 - ratio, 0.);
    const double magnitude =
        outward(outward(2 / (2. * M + 1), infinity) / denominator, infinity);
    const double error = outward(product * magnitude, infinity);
    if (error <= tolerance) return M;
  }
  throw std::runtime_error("truncated Coulomb: zero-seed depth exhausted");
}
}  // namespace

// Budgets contain the two row sums of the contractive propagation branches.
struct QTolerances {
  int join = 0;
  double upward = infinity, budget = 0, downward_h = 0;
  double operator[](unsigned m) const {
    if (m < static_cast<unsigned>(join)) return upward;
    if (m == static_cast<unsigned>(join)) return infinity;
    return budget * ((2. * m + 1) / downward_h);
  }
};

// All mutable storage belongs to the caller and is reused between requests.
struct Workspace {
  std::vector<double> q, origin, hermite, weights, tails;
  std::vector<double> pair_terms;  // origin_pair_taylor
  std::vector<double> inverse_integer, inverse_half;
  QTolerances tolerance;
  std::vector<double> rhs, correction;
  std::vector<std::array<double, 6>> matrix;
  std::vector<unsigned> pivots;
  std::shared_ptr<const FmEval_Chebyshev7<double>> boys;
  int boys_failed_order = INT_MAX;
  double g0 = 0, seed = 0;
  int upward_last = 0, seed_order = -1;
  int origin_last = -1;
  double origin_difference = 0, origin_envelope = 0;
  double exp_minus_t = -1;
  double sqrt_t = -1, displacement = 0;
  bool displacement_ready = false;
  bool gaussian_ready = false;
  double gaussian_minus = 0, gaussian_plus = 0;
  int limit_order = -1, limit_power = 0;
  double limit_precision = 0, limit_gap = 0, limit_min_t = 0;
  void reset_request() {
    origin_last = -1;
    exp_minus_t = -1;
    sqrt_t = -1;
    displacement_ready = false;
    gaussian_ready = false;
    weights.clear();
    tails.clear();
  }
  void ensure_reciprocals(std::size_t size) {
    if (inverse_integer.size() >= size) return;
    static const ExpensiveNumbers<double> numbers;
    const auto begin = inverse_integer.size();
    inverse_integer.resize(size);
    inverse_half.resize(size);
    for (std::size_t k = begin; k < size; ++k) {
      inverse_integer[k] = k ? 1. / k : 0;
      inverse_half[k] = k < 128 ? 2 * numbers.twoi1[k] : 1 / (k + .5);
    }
  }
  explicit Workspace(int n) {
    const std::size_t capacity = static_cast<std::size_t>(n) + 65;
    for (auto* v : {&q, &origin, &hermite, &rhs, &correction, &weights, &tails})
      v->reserve(capacity);
    ensure_reciprocals(capacity);
    matrix.reserve(capacity);
    pivots.reserve(capacity);
    // The standard Engine order range uses exactly the ordinary Coulomb
    // evaluator. Higher radial orders retain an uncapped binary64 limit.
    boys = FmEval_Chebyshev7<double>::instance(0);
    try {
      if (n > boys->max_m()) boys = FmEval_Chebyshev7<double>::instance(n);
    } catch (const std::invalid_argument&) {
      boys_failed_order = n;
    }
  }
};
Scratch::Scratch(int n) {
  if (n < 0)
    throw std::invalid_argument("negative truncated Coulomb scratch order");
  workspace_.reset(new Workspace(n));
}
Scratch::~Scratch() = default;
Scratch::Scratch(const Scratch& other)
    : workspace_(other.workspace_ ? new Workspace(*other.workspace_)
                                  : nullptr) {}
Scratch& Scratch::operator=(const Scratch& other) {
  if (this != &other) {
    if (other.workspace_) {
      if (workspace_)
        *workspace_ = *other.workspace_;
      else
        workspace_.reset(new Workspace(*other.workspace_));
    } else
      workspace_.reset();
  }
  return *this;
}
Scratch::Scratch(Scratch&&) noexcept = default;
Scratch& Scratch::operator=(Scratch&&) noexcept = default;

namespace {
enum class CutoffLimit { none, zero, boys };

// Upper bound on -log F_n(T), hence on -log F_m(T) for all m <= n. For
// T >= n+1/2 the regularized gamma P(n+1/2,T) >= 1/2 (a gamma median is
// below its mean), so F_n >= Gamma(n+1/2)/(4 T^(n+1/2)); otherwise
// F_n >= exp(-T)/(2n+1). The added 1 absorbs rounding of log and lgamma.
double boys_log_scale(double T, int n) {
  const double scale =
      T >= n + .5 ? (n + .5) * std::log(T) - std::lgamma(n + .5) + std::log(4.)
                  : T + std::log(2. * n + 1);
  return std::max(0., scale) + 1;
}

// The gap for an error bound exp(-gap^2) relative to F_n(T).
double relative_gap(double gap, double T, int n) {
  const double square =
      outward_mul_add(gap, gap, boys_log_scale(T, n), infinity);
  return outward(std::sqrt(square), infinity);
}

CutoffLimit finite_cutoff_limit(double R, double T, int n,
                                const Options& options, Workspace& w) {
  // Adjacent binary64 bounds on log(2); no runtime transcendental bound is
  // needed. frexp chooses an exact dyadic error budget <= precision / 32.
  constexpr double log_two_lo = .6931471805599453;
  constexpr double log_two_hi = .6931471805599454;
  if (w.limit_order != n || w.limit_precision != options.precision) {
    int exponent;
    std::frexp(options.precision / 32, &exponent);
    w.limit_power = std::max(0, 1 - exponent);
    const double square =
        outward_mul_add(double(w.limit_power), log_two_hi, 2. * n, infinity);
    w.limit_gap = outward(std::sqrt(square), infinity);
    w.limit_min_t =
        outward_mul_add(double(w.limit_power), log_two_lo, double(n), 0.);
    w.limit_order = n;
    w.limit_precision = options.precision;
  }
  // For a=sqrt(T+n), Cauchy's derivative bound on the disk of radius n
  // gives |F_m-G_m| <= exp(2n-(R-a)^2) when R>=a, simultaneously for
  // m<=n. For n=0 the same inequality follows from the real Gaussian
  // integral. Bound the distance downward and the required gap upward.
  // Both limits are accepted only below precision * F_n(T), not precision
  // alone: at large T the higher orders are far below precision, and an
  // absolute test would replace them by values without a correct digit.
  if (R >= w.limit_gap) {
    const double center = outward(T + n, infinity);
    const double a = outward(std::sqrt(center), infinity);
    if (R >= a) {
      const double distance = outward(R - a, 0.);
      if (distance >= w.limit_gap &&
          distance >= relative_gap(w.limit_gap, T, n))
        return CutoffLimit::boys;
    }
  }
  // In the other direction |G_m| <= R^2 exp(2n-(a-R)^2), R<=a.
  // R^2 <= 1 for R<=1; otherwise frexp gives R^2 < 2^(2*exponent).
  // T >= n + power*log(2) is necessary for this conservative test.
  if (T >= w.limit_min_t) {
    const double center = outward(T + n, 0.);
    const double a = outward(std::sqrt(center), 0.);
    if (a >= R) {
      const double distance = outward(a - R, 0.);
      if (distance >= w.limit_gap) {
        double gap = w.limit_gap;
        if (R > 1) {
          int exponent;
          std::frexp(R, &exponent);
          const double square =
              outward_mul_add(double(w.limit_power + 2 * exponent), log_two_hi,
                              2. * n, infinity);
          gap = outward(std::sqrt(square), infinity);
        }
        if (distance >= gap && distance >= relative_gap(gap, T, n))
          return CutoffLimit::zero;
      }
    }
  }
  return CutoffLimit::none;
}

// Cache coefficients and their tail bounds within one request. G0 and Q
// retain their own first admissible degree, even when tolerances differ.
TaylorPlan taylor_plan(double T, double tolerance, unsigned max_terms,
                       Workspace& w) {
  if (T >= max_terms) return {-1, infinity};
  if (w.weights.empty()) {
    w.weights.push_back(1);
    w.tails.push_back(infinity);
  }
  const double target = tolerance / 128;
  const auto found = std::lower_bound(w.tails.begin() + 1, w.tails.end(),
                                      target, std::greater<double>());
  if (found != w.tails.end())
    return {static_cast<int>(found - w.tails.begin()) - 1, *found};
  if (!std::isfinite(w.weights.back())) return {-1, infinity};
  for (unsigned k = static_cast<unsigned>(w.weights.size()); k <= max_terms;
       ++k) {
    const double weight = w.weights.back() * (-T / k);
    w.weights.push_back(weight);
    const double ratio = T / (k + 1.);
    const double tail = ratio < 1 ? std::abs(weight) / (1 - ratio) : infinity;
    w.tails.push_back(tail);
    if (!std::isfinite(weight)) return {-1, infinity};
    if (tail <= target) return {static_cast<int>(k) - 1, tail};
  }
  return {-1, infinity};
}
// The first omitted G_m coefficient is bounded by
// min(R^2, 2 / (2m + 2k + 1)). Reuse the Q Taylor weights and geometric
// tails, while choosing G's degree from its own, smaller remainder.
TaylorPlan taylor_g_plan(double x, double T, unsigned m, double tolerance,
                         unsigned max_terms, Workspace& w) {
  if (T >= max_terms) return {-1, infinity};
  if (w.weights.empty()) {
    w.weights.push_back(1);
    w.tails.push_back(infinity);
  }
  const double target = tolerance / 128;
  const auto tail_at = [&](std::size_t k) {
    return w.tails[k] * std::min(x, 2 / (2. * m + 2. * k + 1));
  };
  // The cached tails and coefficient bound are both nonincreasing. Find
  // the first admissible degree without rescanning G0's complete prefix.
  std::size_t lo = 1, hi = std::min(w.tails.size(),
                                    static_cast<std::size_t>(max_terms) + 1);
  const std::size_t end = hi;
  while (lo < hi) {
    const auto mid = lo + (hi - lo) / 2;
    if (tail_at(mid) <= target)
      hi = mid;
    else
      lo = mid + 1;
  }
  if (lo < end) return {static_cast<int>(lo) - 1, tail_at(lo)};
  if (!std::isfinite(w.weights.back())) return {-1, infinity};
  for (unsigned k = static_cast<unsigned>(w.weights.size()); k <= max_terms;
       ++k) {
    const double weight = w.weights.back() * (-T / k);
    w.weights.push_back(weight);
    const double ratio = T / (k + 1.);
    w.tails.push_back(ratio < 1 ? std::abs(weight) / (1 - ratio) : infinity);
    if (!std::isfinite(weight)) return {-1, infinity};
    const double tail = tail_at(k);
    if (tail <= target) return {static_cast<int>(k) - 1, tail};
  }
  return {-1, infinity};
}
double exponential(double T, Workspace& w) {
  if (w.exp_minus_t < 0) w.exp_minus_t = std::exp(-T);
  return w.exp_minus_t;
}
double square_root(double T, Workspace& w) {
  if (w.sqrt_t < 0) w.sqrt_t = std::sqrt(T);
  return w.sqrt_t;
}
double displacement(double R, double T, Workspace& w) {
  if (!w.displacement_ready) {
    w.displacement = root_difference(R, T, square_root(T, w));
    w.displacement_ready = true;
  }
  return w.displacement;
}
void gaussians(double R, double T, Workspace& w) {
  if (w.gaussian_ready) return;
  const double z = displacement(R, T, w), plus = R + square_root(T, w);
  w.gaussian_minus = std::exp(-z * z);
  w.gaussian_plus = std::exp(-plus * plus);
  w.gaussian_ready = true;
}
void boys_limit(double T, int n, Workspace& w, double* g) {
  // A workspace reserved beyond the native table can later receive a lower
  // request within it. Select by the actual order, independently of reserve
  // capacity. Remember failed growth so unsupported orders do not keep
  // throwing.
  if (n > w.boys->max_m() && n < w.boys_failed_order) {
    try {
      w.boys = FmEval_Chebyshev7<double>::instance(n);
    } catch (const std::invalid_argument&) {
      w.boys_failed_order = n;
    }
  }
  if (w.boys && n <= w.boys->max_m()) {
    w.boys->eval(g, T, n);
  } else {
    // Beyond the existing Boys table, join two contractive recurrences.
    // In particular, an underflowed exp(-T) must not erase F0, and an
    // upward recurrence must not grow a parasitic solution above its turn.
    w.boys->eval(g, T, 0);
    const double e = std::exp(-T);
    int upward = 0;
    while (upward < n && upward + .5 < T) {
      g[upward + 1] = mul_add(2. * upward + 1, g[upward], -e) / (2 * T);
      ++upward;
    }
    if (upward < n) {
      double value = 0;
      if (e != 0) {
        double term = 1 / (2. * n + 1);
        Sum sum;
        sum.add(term);
        for (int k = 1;; ++k) {
          term *= T / (n + k + .5);
          sum.add(term);
          if (term <= unit * sum.value) break;
        }
        value = e * sum.value;
      }
      g[n] = value;
      for (int m = n - 1; m > upward; --m) {
        value = mul_add(2 * T, value, e) / (2. * m + 1);
        g[m] = value;
      }
    }
  }
}
void evaluate_limit(double T, int n, CutoffLimit limit, Workspace& w,
                    double* g) {
  if (limit == CutoffLimit::boys)
    boys_limit(T, n, w, g);
  else
    std::fill(g, g + static_cast<std::size_t>(n) + 1, 0);
}
void origin_q(double R, unsigned last, Workspace& w) {
  if (w.origin_last >= static_cast<int>(last)) return;
  w.origin.resize(static_cast<std::size_t>(last) + 1);
  const double x = R * R;
  if (!std::isfinite(x)) {
    std::fill(w.origin.begin(), w.origin.end(), 0);
    w.origin_envelope = 0;
    w.origin_last = static_cast<int>(last);
    return;
  }
  if (w.origin_last < 0) {
    w.origin[0] = std::exp(-x);
    w.origin_envelope = w.origin[0];
    w.origin_difference = 0;
    w.origin_last = 0;
  }
  w.ensure_reciprocals(static_cast<std::size_t>(last));
  for (unsigned k = static_cast<unsigned>(w.origin_last); k < last; ++k) {
    w.origin_difference =
        mul_add(k, w.origin_difference, -x * w.origin[k]) * w.inverse_half[k];
    w.origin[k + 1] = w.origin[k] + w.origin_difference;
    w.origin_envelope = std::max(w.origin_envelope, std::abs(w.origin[k + 1]));
  }
  w.origin_last = static_cast<int>(last);
}
Candidate r_series(double R, double T, unsigned m, bool g, double tolerance,
                   unsigned max_terms, Workspace& w) {
  const double x = R * R, e = exponential(T, w);
  if (!std::isfinite(x) || e == 0) return {0, infinity};
  double previous = 0, term = g ? x * e : e, majorant = std::abs(term);
  const double alpha1 = g ? m + 1.5 : m + .5;
  double ratio = x * (alpha1 + T) / (g ? 3. : .5);
  Sum sum;
  for (unsigned j = 0; j < max_terms; ++j) {
    sum.add(term);
    const double next_majorant = majorant * ratio;
    const double next_denominator =
        g ? (j + 3.) * (j + 2.5) : (j + 2.) * (j + 1.5);
    const double next_ratio = x * (2 * j + 2 + alpha1 + T) / next_denominator;
    if (!std::isfinite(next_majorant)) return {0, infinity};
    if (next_ratio < 1) {
      const double tail = next_majorant / (1 - next_ratio);
      if (tail <= tolerance / 128)
        return {sum.value, 16 * unit * (j + 1) * sum.magnitude + tail};
    }
    // No signed next term is needed after convergence or resource exhaustion.
    if (j + 1 == max_terms) break;
    const double denominator = g ? (j + 2.) * (j + 1.5) : (j + 1.) * (j + .5);
    const double a = -x * (2 * j + alpha1 - T) / denominator;
    const double b = j == 0 ? 0
                     : g    ? (m == 0 ? x * x * j / (denominator * (j + 1.))
                                      : x * x * j * (j + m + .5) /
                                         (denominator * (j + 1.) * (j + .5)))
                            : x * x * (j + m - .5) / (denominator * (j - .5));
    const double next = mul_add(a, term, -b * previous);
    if (!std::isfinite(next)) return {0, infinity};
    previous = term;
    term = next;
    majorant = next_majorant;
    ratio = next_ratio;
  }
  return {0, infinity};
}
Candidate taylor_g0(double R, const TaylorPlan& plan, Workspace& w) {
  origin_q(R, static_cast<unsigned>(plan.degree), w);
  Sum sum;
  sum.add(-std::expm1(-R * R));
  for (int j = 1; j <= plan.degree; ++j)
    sum.add(w.weights[j] * (1 - w.origin[j]) / (2. * j + 1));
  return {sum.value, 8 * unit * sum.magnitude + plan.remainder};
}
Candidate direct_seed(double R, double T, unsigned m, const Options& options,
                      Workspace& w) {
  const double x = R * R;
  // The same conditioning bound as G0's radius series, at the seed order.
  if (x * ((m + 1.5 + T) / 3) <= .25) {
    const auto value =
        r_series(R, T, m, true, options.precision, options.max_taylor_terms, w);
    if (std::isfinite(value.error)) return value;
  }
  // The origin coefficients obey 0 <= G_k(R,0) <= 2/(2k+1).
  // exp(T) <= 2 bounds the Taylor absolute-term sum, as for G0.
  if (T <= std::log(2.)) {
    const auto plan =
        taylor_g_plan(x, T, m, options.precision, options.max_taylor_terms, w);
    if (plan.degree >= 0) {
      origin_q(R, m + static_cast<unsigned>(plan.degree), w);
      Sum sum;
      for (int j = 0; j <= plan.degree; ++j)
        sum.add(w.weights[j] * (1 - w.origin[m + j]) / (2. * m + 2. * j + 1));
      if (std::isfinite(sum.value))
        return {sum.value, 8 * unit * sum.magnitude + plan.remainder};
    }
  }
  return {0, infinity};
}
Candidate evaluate_g0(double R, double T, const Options& options,
                      Workspace& w) {
  if (T == 0) {
    const double value = -std::expm1(-R * R);
    return {value, 4 * unit * std::abs(value)};
  }
  const double x = R * R;
  // r_0 <= 1/4 bounds the positive majorant / |sum| by two.
  if (x * ((T + 1.5) / 3) <= .25) {
    const auto value =
        r_series(R, T, 0, true, options.precision, options.max_taylor_terms, w);
    if (std::isfinite(value.error)) return value;
  }
  // For exp(T) <= 2 the origin Taylor series has a bounded absolute-term
  // sum. It also avoids dividing rounded erf(sqrt(T)) by a tiny sqrt(T).
  if (T <= std::log(2.)) {
    const auto plan =
        taylor_plan(T, options.precision, options.max_taylor_terms, w);
    if (plan.degree >= 0) return taylor_g0(R, plan, w);
  }
  const double s = square_root(T, w), z = displacement(R, T, w);
  const double factor = .44311346272637900682 / s;
  double value, magnitude;
  if (R >= s) {
    double f;
    w.boys->eval(&f, T, 0);
    const double a = factor * std::erfc(-z), b = factor * std::erfc(R + s);
    value = f - a + b;
    magnitude = f + a + b;
    // Near R = s this form cancels by up to 3. With F0 = 2 factor erf(s)
    // the difference form below is the same function, and it is well
    // conditioned once erfc(s) is small.
    if (magnitude > 2 * std::abs(value)) {
      const double c = std::erfc(z), d = std::erfc(s + R), e = 2 * std::erfc(s);
      if (factor * (c + d + e) < magnitude) {
        value = factor * (c + d - e);
        magnitude = factor * (c + d + e);
      }
    }
  } else {
    const double a = std::erfc(z), b = std::erfc(s + R), c = 2 * std::erfc(s);
    value = factor * (a + b - c);
    magnitude = factor * (a + b + c);
  }
  if (magnitude <= 2 * std::abs(value) || magnitude <= options.precision / 16)
    return {value, 16 * unit * magnitude};
  const auto radius =
      r_series(R, T, 0, true, options.precision, options.max_taylor_terms, w);
  if (std::isfinite(radius.error) && radius.error < 16 * unit * magnitude)
    return radius;
  const auto plan =
      taylor_plan(T, options.precision, options.max_taylor_terms, w);
  if (plan.degree >= 0) {
    const auto taylor = taylor_g0(R, plan, w);
    if (std::isfinite(taylor.error) && taylor.error < 16 * unit * magnitude)
      return taylor;
  }
  return {value, 16 * unit * magnitude};
}

void gaussian_q0(double R, double T, Workspace& w) {
  gaussians(R, T, w);
  // This is also Miller's known boundary, so solving cannot improve it.
  w.q[0] = .5 * (w.gaussian_minus + w.gaussian_plus);
}
void scaled_hermites(double R, double T, unsigned last, Workspace& w) {
  auto& p = w.hermite;
  p.resize(static_cast<std::size_t>(last) + 1);
  gaussians(R, T, w);
  const bool minus_active = w.gaussian_minus != 0;
  const bool plus_active = w.gaussian_plus != 0;
  p[0] = w.gaussian_minus + w.gaussian_plus;
  if (!minus_active && !plus_active) {
    std::fill(p.begin(), p.end(), 0);
    return;
  }
  if (!last) return;
  const double s = square_root(T, w);
  const double a_minus = minus_active ? displacement(R, T, w) / s : 0;
  const double a_plus = plus_active ? (s + R) / s : 0;
  double minus_previous = w.gaussian_minus, plus_previous = w.gaussian_plus;
  double minus = a_minus * minus_previous, plus = a_plus * plus_previous;
  p[1] = minus + plus;
  const double inv2T = .5 / T;
  for (unsigned k = 1; k < last; ++k) {
    const double b = k * inv2T;
    const double next_minus =
        minus_active ? mul_add(a_minus, minus, -b * minus_previous) : 0;
    const double next_plus =
        plus_active ? mul_add(a_plus, plus, -b * plus_previous) : 0;
    p[k + 1] = next_minus + next_plus;
    minus_previous = minus;
    plus_previous = plus;
    minus = next_minus;
    plus = next_plus;
  }
}

bool hermite_q(double R, double T, unsigned last, Workspace& w) {
  scaled_hermites(R, T, last, w);
  w.ensure_reciprocals(static_cast<std::size_t>(last) + 1);
  // Dividing .25 by T avoids overflowing a denominator 4*T at large T.
  const double inv4T = .25 / T;
  for (unsigned m = 1; m <= last; ++m) {
    if (m == static_cast<unsigned>(w.upward_last)) continue;
    Sum sum;
    double weight = .5;
    for (unsigned k = m; k > 0; --k) {
      sum.add(weight * w.hermite[k]);
      if (k > 1)
        weight *=
            (2. * m - k) * (k - 1.) * inv4T * w.inverse_integer[m - k + 1];
    }
    if (!std::isfinite(sum.value) || std::abs(sum.value) > 1 ||
        unit * std::max(0., sum.magnitude - std::abs(sum.value)) >
            w.tolerance[m])
      return false;
    w.q[m] = sum.value;
  }
  return true;
}

double taylor_absolute_weights(const TaylorPlan& plan, double gamma,
                               const Workspace& w) {
  if (plan.degree == 0) return 1;
  const double normal = std::numeric_limits<double>::min();
  // |w_j| is unimodal. These endpoint tests ensure normal intermediates
  // in the two-operation coefficient recurrence for every j <= degree.
  if (std::abs(w.weights[plan.degree]) >= normal &&
      std::abs(w.weights[1]) >= normal * plan.degree && w.exp_minus_t > 0) {
    if (taylor_roundings == 1) {
      const double du = (2. * plan.degree) * unit;
      gamma = outward(du / (1 - du), infinity);
    }
    const double denominator = outward(w.exp_minus_t, 0.);
    if (denominator > 0)
      return outward(outward(1 + gamma, infinity) / denominator, infinity);
  }
  // Subnormal coefficients need an absolute bound on the stored prefix.
  double result = 1;
  for (int j = 1; j <= plan.degree; ++j)
    result = outward(result + std::abs(w.weights[j]), infinity);
  return result;
}

bool taylor_q(double R, unsigned last, const TaylorPlan& plan, Workspace& w) {
  origin_q(R, last + static_cast<unsigned>(plan.degree), w);
  const unsigned low_end =
      std::min(last + 1, static_cast<unsigned>(w.upward_last));
  const unsigned high_begin =
      std::min(last + 1, static_cast<unsigned>(w.upward_last) + 1);
  unsigned fast_low_end = 1, fast_high_begin = last + 1;
  const double operations = taylor_roundings * static_cast<double>(plan.degree);
  const double du = operations * unit;
  const double gamma = du == 0 ? 0 : outward(du / (1 - du), infinity);
  double largest_tolerance = low_end > 1 ? w.tolerance[1] : 0;
  if (high_begin <= last)
    largest_tolerance = std::max(largest_tolerance, w.tolerance[last]);
  // A cheap exp(-T) predictor avoids constructing an accumulation bound
  // when no requested order could use ordinary summation. The final test
  // below uses the actual stored weights and origin envelope, rather than
  // treating rounded coefficients as exact mathematical coefficients.
  if (gamma * w.origin_envelope <=
      (largest_tolerance - plan.remainder) * w.exp_minus_t) {
    const double absolute_weights = taylor_absolute_weights(plan, gamma, w);
    const double magnitude_bound =
        outward(absolute_weights * w.origin_envelope, infinity);
    const double underflow = outward(
        operations * std::numeric_limits<double>::denorm_min() / (1 - du),
        infinity);
    const double accumulation_error = outward(
        outward(gamma * magnitude_bound, infinity) + underflow, infinity);
    const double error = outward(accumulation_error + plan.remainder, infinity);
    if (low_end > 1 && error <= w.tolerance[1]) fast_low_end = low_end;
    // Downward tolerances increase with m, so the eligible orders form a
    // single contiguous suffix. The unused join separates the two ranges.
    unsigned lo = high_begin, hi = last + 1;
    while (lo < hi) {
      const unsigned mid = lo + (hi - lo) / 2;
      if (error <= w.tolerance[mid])
        hi = mid;
      else
        lo = mid + 1;
    }
    fast_high_begin = lo;
  }
  std::copy(w.origin.begin() + 1, w.origin.begin() + fast_low_end,
            w.q.begin() + 1);
  std::copy(w.origin.begin() + fast_high_begin, w.origin.begin() + last + 1,
            w.q.begin() + fast_high_begin);
  // Independent m lanes permit vector multiply/add, or native vector FMA.
  if (fast_low_end > 1 || fast_high_begin <= last) {
    for (int j = 1; j <= plan.degree; ++j) {
      const double weight = w.weights[j];
      for (unsigned m = 1; m < fast_low_end; ++m)
        w.q[m] = taylor_add(weight, w.origin[m + j], w.q[m]);
      for (unsigned m = fast_high_begin; m <= last; ++m)
        w.q[m] = taylor_add(weight, w.origin[m + j], w.q[m]);
    }
  }
  for (unsigned m = 1; m <= last; ++m) {
    if (m == static_cast<unsigned>(w.upward_last)) continue;
    if (m < fast_low_end || m >= fast_high_begin) {
      if (!std::isfinite(w.q[m]) || std::abs(w.q[m]) > 1 + unit * 2)
        return false;
      continue;
    }
    Sum sum;
    sum.add(w.origin[m]);
    for (int j = 1; j <= plan.degree; ++j) {
      sum.add(w.weights[j] * w.origin[m + j]);
    }
    if (!std::isfinite(sum.value) || std::abs(sum.value) > 1 + unit * 2)
      return false;
    w.q[m] = sum.value;
  }
  return true;
}
// Extending the artificial boundary reuses the outward recurrence prefix.
// At an accepted endpoint its product has already included that order,
// while the interval state still precedes the next recurrence step.
struct MillerBoundary {
  double T, xlo, xhi;
  unsigned last, order = 0;
  double qlo = 0, qhi = 0, hlo = 0, hhi = 0, product = 1;
  bool counted = false, failed = false;
  MillerBoundary(double R, double t, unsigned l) : T(t), last(l) {
    const double x = R * R;
    xlo = outward(x, 0.);
    xhi = outward(x, infinity);
  }
  double advance(unsigned n) {
    if (failed) return infinity;
    for (;;) {
      if (!counted && order >= last) {
        const double denominator = outward(qhi + 1, 0.);
        const double ratio = std::min(1., outward(qhi / denominator, infinity));
        product = outward(product * ratio, infinity);
      }
      counted = true;
      if (order == n) return product;
      const double flo = outward_mul_add(double(order), hlo,
                                         outward_mul_add(xlo, qlo, .5, 0.), 0.);
      const double fhi =
          outward_mul_add(double(order), hhi,
                          outward_mul_add(xhi, qhi, .5, infinity), infinity);
      const double dlo = outward(T + flo, 0.), dhi = outward(T + fhi, infinity);
      qlo = outward(outward(T * outward(qlo + 1, 0.), 0.) / dhi, 0.);
      qhi = outward(outward(T * outward(qhi + 1, infinity), infinity) / dlo,
                    infinity);
      // Retain the shared f in the monotone ratio f/(T+f).
      hlo = outward(flo / outward(T + flo, infinity), 0.);
      hhi = outward(fhi / outward(T + fhi, 0.), infinity);
      if (!std::isfinite(qhi) || !std::isfinite(hhi)) {
        failed = true;
        return infinity;
      }
      ++order;
      counted = false;
    }
  }
};

bool miller_q(double R, double T, unsigned last, double tolerance,
              unsigned max_extra, Workspace& w) {
  if (T == 0 || !std::isfinite(R * R)) return false;
  gaussians(R, T, w);
  const double q0 = .5 * (w.gaussian_minus + w.gaussian_plus);
  if (q0 == 0) return false;
  max_extra = std::min(max_extra, static_cast<unsigned>(INT_MAX) - last);
  if (!max_extra) return false;
  MillerBoundary boundary_state(R, T, last);
  for (unsigned extra = std::min(16u, max_extra); extra > 0;) {
    const unsigned n = last + extra;
    const double boundary = boundary_state.advance(n);
    if (!(boundary <= tolerance / 128)) {
      if (extra == max_extra) break;
      extra += std::min(extra, max_extra - extra);
      continue;
    }
    // Column j of row i is stored at j-i+2. Fill-in after pivoting extends
    // through column i+3; no dense matrix or external linear algebra is used.
    auto& matrix = w.matrix;
    matrix.resize(n);
    auto& rhs = w.rhs;
    rhs.resize(n);
    auto& pivots = w.pivots;
    pivots.resize(n);
    for (unsigned m = 0; m < n; ++m) {
      matrix[m].fill(0);
      rhs[m] = 0;
      const double coefficients[] = {-double(m), T + 2 * m + .5 - R * R,
                                     -(2 * T + m + .5), T};
      for (unsigned k = 0; k < 4; ++k) {
        const int order = static_cast<int>(m) + static_cast<int>(k) - 1;
        if (order == 0) rhs[m] -= coefficients[k];
        if (order > 0 && order <= static_cast<int>(n))
          matrix[m][k] = coefficients[k];
      }
    }
    bool finite = true;
    for (unsigned k = 0; k < n; ++k) {
      unsigned pivot = k;
      for (unsigned i = k + 1; i < std::min(n, k + 3); ++i)
        if (std::abs(matrix[i][k - i + 2]) >
            std::abs(matrix[pivot][k - pivot + 2]))
          pivot = i;
      if (pivot != k) {
        for (unsigned j = k; j < std::min(n, k + 4); ++j)
          std::swap(matrix[k][j - k + 2], matrix[pivot][j - pivot + 2]);
        std::swap(rhs[k], rhs[pivot]);
      }
      pivots[k] = pivot;
      const double diagonal = matrix[k][2];
      if (diagonal == 0 || !std::isfinite(diagonal)) {
        finite = false;
        break;
      }
      for (unsigned i = k + 1; i < std::min(n, k + 3); ++i) {
        const double factor = matrix[i][k - i + 2] / diagonal;
        for (unsigned j = k + 1; j < std::min(n, k + 4); ++j)
          matrix[i][j - i + 2] =
              mul_add(-factor, matrix[k][j - k + 2], matrix[i][j - i + 2]);
        rhs[i] = mul_add(-factor, rhs[k], rhs[i]);
        matrix[i][k - i + 2] = factor;
      }
    }
    if (!finite) return false;
    for (unsigned kk = n; kk > 0; --kk) {
      const unsigned k = kk - 1;
      for (unsigned j = k + 1; j < std::min(n, k + 4); ++j)
        rhs[k] = mul_add(-matrix[k][j - k + 2], rhs[j], rhs[k]);
      rhs[k] /= matrix[k][2];
    }
    for (auto& value : rhs) value *= q0;
    auto& correction = w.correction;
    correction.resize(n);
    for (unsigned iteration = 0; iteration < 2; ++iteration) {
      // Form the residual in differences, before collecting the diagonal
      // coefficient. This retains R^2 when 2m is much larger and preserves
      // a nearly constant hierarchy. Refinement uses the same double LU.
      for (unsigned m = 0; m < n; ++m) {
        const double current = m == 0 ? q0 : rhs[m - 1];
        const double previous = m == 0 ? 0 : m == 1 ? q0 : rhs[m - 2];
        const double next = rhs[m], after = m + 1 == n ? 0 : rhs[m + 1];
        const double difference = next - current;
        const double residual = mul_add(
            T, (after - next) - difference,
            mul_add(double(m), current - previous, -(m + .5) * difference));
        correction[m] = mul_add(R * R, current, -residual);
      }
      for (unsigned k = 0; k < n; ++k) {
        std::swap(correction[k], correction[pivots[k]]);
        for (unsigned i = k + 1; i < std::min(n, k + 3); ++i)
          correction[i] =
              mul_add(-matrix[i][k - i + 2], correction[k], correction[i]);
      }
      for (unsigned kk = n; kk > 0; --kk) {
        const unsigned k = kk - 1;
        for (unsigned j = k + 1; j < std::min(n, k + 4); ++j)
          correction[k] =
              mul_add(-matrix[k][j - k + 2], correction[j], correction[k]);
        correction[k] /= matrix[k][2];
      }
      for (unsigned k = 0; k < n; ++k) rhs[k] += correction[k];
    }
    double envelope = q0;
    for (unsigned m = 1; m <= last; ++m)
      envelope = std::max(envelope, std::abs(rhs[m - 1]));
    if (std::isfinite(envelope)) {
      w.q[0] = q0;
      for (unsigned m = 1; m <= last; ++m) w.q[m] = rhs[m - 1];
      return true;
    }
    if (extra == max_extra) break;
    extra += std::min(extra, max_extra - extra);
  }
  return false;
}

// Uniform injection budgets, scaled by the exact recurrence coefficients.
// Along each contractive branch the coefficients increase in traversal
// order, so its positive row sums increase up to the join. H is the final
// row sum; no per-order tolerance array or maximum reduction is needed.
// Thus sum_j |d G_m / d Q_j| tolerance_j <= budget for every requested m,
// including a join arbitrarily close to an integer-plus-one-half.
double q_tolerances(double T, int last, double budget, Workspace& w) {
  auto& tolerance = w.tolerance;
  tolerance.join = w.upward_last;
  tolerance.budget = budget;
  double h = 0;
  for (int m = 0; m < w.upward_last; ++m) {
    h = mul_add((m + .5) / T, h, 1.);
  }
  tolerance.upward = h > 0 ? budget * (2 * T / h) : infinity;
  h = 0;
  for (int m = last; m > w.upward_last; --m) {
    h = mul_add(T / (m + .5), h, 1.);
  }
  tolerance.downward_h = h;
  double minimum = tolerance.join > 1 ? tolerance.upward : infinity;
  if (last > tolerance.join)
    minimum = std::min(minimum, tolerance[tolerance.join + 1]);
  return minimum;
}

void evaluate_q(double R, double T, unsigned last, double tolerance,
                const Options& options, Workspace& w) {
  w.q.resize(static_cast<std::size_t>(last) + 1);
  if (T == 0) {
    origin_q(R, last, w);
    std::copy(w.origin.begin(), w.origin.begin() + last + 1, w.q.begin());
    return;
  }
  // An unused Q0 needs no shifted Gaussians on either series path.
  if (w.upward_last > 0) gaussian_q0(R, T, w);
  if (last == 0) return;
  const double x = R * R;
  if (x * (2 * (last + .5 + T)) <= .25) {
    bool finite = true;
    for (unsigned m = 1; m <= last; ++m) {
      if (m == static_cast<unsigned>(w.upward_last)) continue;
      const double target =
          std::isfinite(w.tolerance[m]) ? w.tolerance[m] : options.precision;
      const auto value =
          r_series(R, T, m, false, target, options.max_taylor_terms, w);
      if (!std::isfinite(value.error)) {
        finite = false;
        break;
      }
      w.q[m] = value.value;
    }
    if (finite) return;
  }
  // |Q_m(R,0)| <= 1 bounds the Taylor absolute-term sum by exp(T).
  // This predicts conditioning without first evaluating every candidate.
  TaylorPlan plan{-1, infinity};
  bool taylor_tried = false;
  if (unit <= tolerance * exponential(T, w)) {
    taylor_tried = true;
    plan = taylor_plan(T, tolerance, options.max_taylor_terms, w);
    if (plan.degree >= 0 && taylor_q(R, last, plan, w)) return;
  }
  if (hermite_q(R, T, last, w)) return;
  if (miller_q(R, T, last, tolerance, options.max_extra_orders, w)) return;
  // Resource exhaustion or underflow in one representation is a fallback,
  // rather than a reason to enumerate all representations on every call.
  if (!taylor_tried) {
    plan = taylor_plan(T, tolerance, options.max_taylor_terms, w);
    if (plan.degree >= 0 && taylor_q(R, last, plan, w)) return;
  }
  throw std::runtime_error(
      "truncated Coulomb: no finite bounded Q representation");
}

// Double-word tables of 1/((1/2)_k (2k+1)) and 1/k!. Each entry takes at
// most k double-word scalings by a double, each with a relative error of a
// few u^2 (Joldes, Muller, Popescu, ACM TOMS 44, 2017), and one Newton
// reciprocal. The bound 512 u^2 = 2^-97 used below exceeds the 9 u^2
// measured against 60-digit values; FMA and Dekker products give the same
// entries.
constexpr int origin_pair_capacity = 160;
constexpr double origin_pair_max_t = 12;
struct OriginPairTables {
  std::array<Pair, origin_pair_capacity> coefficient, inverse_factorial;
  // log (1/2)_k, and kappa_k = k! / (1/2)_{k+1} for the error estimate
  std::array<double, origin_pair_capacity> log_pochhammer, kappa;
  static constexpr double error = 6.310887241768094e-30;  // 2^-97
  OriginPairTables() {
    Pair pochhammer{1, 0}, factorial{1, 0};
    for (int k = 0; k < origin_pair_capacity; ++k) {
      coefficient[k] = reciprocal(scale(pochhammer, 2. * k + 1));
      inverse_factorial[k] = reciprocal(factorial);
      log_pochhammer[k] = std::log(pochhammer.hi);
      pochhammer = scale(pochhammer, k + .5);
      kappa[k] = factorial.hi / pochhammer.hi;
      factorial = scale(factorial, k + 1.);
    }
  }
  static Pair scale(Pair a, double b) {
    const Pair p = two_product(a.hi, b);
    return fast_two_sum(p.hi, p.lo + a.lo * b);
  }
  static Pair reciprocal(Pair d) {
    const double r = 1 / d.hi;
    const Pair p = two_product(d.hi, r);
    return fast_two_sum(r, r * (((1 - p.hi) - p.lo) - d.lo * r));
  }
};

// Taylor expansion about T = 0 for all orders at once, in double-word
// arithmetic: with x = R^2,
//   G_m = F_m(T) - exp(-x) S_m,  S_m = sum_j (-T)^j / j! p_{m+j},
// p_k = P_k(x) / (2k+1) and P_k = 1F1(-k; 1/2; x). Q_k(R,0) = exp(-x) P_k,
// so |exp(-x) p_k| <= 1/(2k+1) (see direct_seed). The absolute terms of
// exp(-x) S_m sum to at most exp(T)/(2m+1), so binary64 summation would lose
// log10 exp(T) digits; the double-word error stays near u^2 exp(T).
// U_k = (1/2)_k P_k obeys the division-free recurrence
//   U_{k+1} = c_k U_k + d_k U_{k-1},  c_k = 2k + 1/2 - x,  d_k = -k(k - 1/2),
// and p_k = U_k / ((1/2)_k (2k+1)) uses a table. Returns false, leaving g
// unspecified, when the error estimate exceeds the budget or the request
// exceeds the tables.
bool origin_pair_taylor(double R, double T, int n, const Options& options,
                        Workspace& w, double* g) {
  static const OriginPairTables tables;
  constexpr int capacity = origin_pair_capacity;
  // sum_{j>d} T^j / (j! (2j+1)) bounds the truncation of exp(-x) S_m.
  const double target = options.precision / 128;
  double term = 1, tail = infinity;
  int degree = 0;
  for (; n + degree + 1 < capacity; ++degree) {
    term *= T / (degree + 1);
    const double ratio = T / (degree + 2.);
    if (ratio < 1) {
      tail = term / ((1 - ratio) * (2. * degree + 3));
      if (tail <= target) break;
    }
  }
  if (!(tail <= target)) return false;
  const int last = n + degree;
  const Pair x = two_product(R, R);
  // |U_k| <= (1/2)_k exp(x) < 2^900 keeps every split exact, and the
  // error estimate below finite.
  if (!(x.hi <= 300) || tables.log_pochhammer[last] + x.hi > 600) return false;
  const double ex = std::exp(-x.hi) * (1 - x.lo);
  const double underflow = 16 * std::numeric_limits<double>::denorm_min();

  w.pair_terms.resize(14 * capacity);
  double* const p_hi = w.pair_terms.data();
  double* const p_lo = p_hi + capacity;
  double* const p_split_hi = p_lo + capacity;
  double* const p_split_lo = p_split_hi + capacity;
  double* const w_hi = p_split_lo + capacity;
  double* const w_lo = w_hi + capacity;
  double* const s_hi = w_lo + capacity;
  double* const s_lo = s_hi + capacity;
  double* const u_hi = s_lo + capacity;  // U_k
  double* const u_lo = u_hi + capacity;
  double* const v = u_lo + capacity;  // V_k
  double* const alpha = v + capacity;
  double* const P = alpha + capacity;  // |exp(-x) P_k|
  double* const Q = P + capacity;      // |Q_k|

  // U_k, with c_k = c.hi + c_lo (the latter rounded once) and exact d_k.
  // tau_{k+1} bounds the rounding in forming U_{k+1}, to first order: the
  // products c.hi U_k.lo and c_lo U_k.hi and the sums after them (5 times),
  // d U_{k-1}.lo (4 times), the error-free product and sum residuals (u
  // times their leading parts, 3 and 1 times), the rounding of c_lo and the
  // dropped c_lo U_k.lo. alpha_{k+1} = tau_{k+1} / k! is used below; it is
  // independent of the recurrence and costs little here.
  const double u = unit;
  const Pair c0 = two_sum(.5, -x.hi);
  u_hi[0] = 1;
  u_lo[0] = 0;
  u_hi[1] = c0.hi;
  u_lo[1] = c0.lo - x.lo;
  v[0] = 0;
  v[1] = 1;
  alpha[1] = u * std::abs(u_lo[1]) + underflow;
  double magnitude_previous = 3 * u;  // 5 |U_0.lo| + 3 u |U_0.hi|
  for (int k = 1; k < last; ++k) {
    const Pair c = two_sum(2. * k + .5, -x.hi);
    const double c_lo = c.lo - x.lo, d = -k * (k - .5);
    const Pair p1 = two_product(c.hi, u_hi[k]);
    const Pair p2 = two_product(d, u_hi[k - 1]);
    const double e1 = p1.lo + (c.hi * u_lo[k] + c_lo * u_hi[k]);
    const double e2 = p2.lo + d * u_lo[k - 1];
    const Pair sum = two_sum(p1.hi, p2.hi);
    u_hi[k + 1] = sum.hi;
    u_lo[k + 1] = sum.lo + (e1 + e2);
    v[k + 1] = mul_add(c.hi, v[k], d * v[k - 1]);
    const double magnitude = 5 * std::abs(u_lo[k]) + 3 * u * std::abs(u_hi[k]);
    const double tau =
        u * (std::abs(c.hi) * magnitude - d * magnitude_previous +
             u * std::abs(sum.hi)) +
        std::abs(c_lo) * (std::abs(u_lo[k]) + 6 * u * std::abs(u_hi[k])) +
        underflow;
    alpha[k + 1] = tau * tables.inverse_factorial[k].hi;
    magnitude_previous = magnitude;
  }
  // p_k = U_k c_k; P_k and Q_k = V_k / (1/2)_k for the estimate below.
  p_hi[0] = 1;
  p_lo[0] = 0;
  P[0] = ex;
  Q[0] = 0;
  double max_p = 1, max_p_lo = 0;
  for (int k = 1; k <= last; ++k) {
    const Pair coefficient = tables.coefficient[k];
    const Pair product = two_product(u_hi[k], coefficient.hi);
    p_hi[k] = product.hi;
    p_lo[k] =
        product.lo + (u_hi[k] * coefficient.lo + u_lo[k] * coefficient.hi);
    const double order = 2. * k + 1;
    P[k] = ex * order * std::abs(product.hi);
    Q[k] = order * coefficient.hi * std::abs(v[k]);
    max_p = std::max(max_p, std::abs(product.hi));
    max_p_lo = std::max(max_p_lo, std::abs(p_lo[k]));
  }
  // Rounding of p_k: with |c_k.lo| <= u c_k.hi the three roundings after
  // the exact product are at most u (3 |p_k.lo| + 10 u |p_k.hi|), to first
  // order; the table adds its relative error.
  const double p_rounding =
      ex * (u * (3 * max_p_lo + 10 * u * max_p) + tables.error * max_p) +
      underflow;

  // An estimate of the error of exp(-x) P_k = (2k+1) exp(-x) p_k. The
  // recurrence error follows its Green's function: a local error tau_i in
  // U_i reaches U_k as
  //   tau_i (U_{i-1} V_k - V_{i-1} U_k) / C_{i-1},
  // with V a second solution (V_0 = 0, V_1 = 1) and the Casoratian
  // C_i = U_i V_{i+1} - U_{i+1} V_i = i! (1/2)_i. With kappa_i =
  // i! / (1/2)_{i+1} the error of exp(-x) P_k is at most
  // sum_i alpha_i |g_ki|, where
  //   g_ki = exp(-x) (P_{i-1} Q_k - Q_{i-1} P_k)
  //        = exp(-x) Q_{i-1} Q_k sum_{l=i-1}^{k-1} kappa_l / (Q_l Q_{l+1}).
  // Bounding the first form term by term suits oscillating U and V. Where
  // both grow with the dominant solution (x >~ 4k) its terms cancel; the
  // second form, as prefix sums of positive terms, then stays near the
  // true error, and is evaluated when the first one misses the budget.
  // This is a first-order estimate: V and the sums are rounded binary64,
  // and the result enters doubled below.
  const auto recurrence_error = [&](bool nested_form) -> double {
    double sum_P = 0, sum_Q = 0, nested = 0, result = 0;
    const double first = ex * alpha[1];  // the i = 1 term; Q_0 = 0
    for (int k = 1; k <= last; ++k) {
      sum_P += P[k - 1] * alpha[k];
      sum_Q += Q[k - 1] * alpha[k];
      double bound = Q[k] * sum_P + P[k] * sum_Q;
      if (nested_form) {
        if (k > 1) nested += tables.kappa[k - 1] / (Q[k - 1] * Q[k]) * sum_Q;
        bound = std::min(bound, Q[k] * (first + ex * nested));
      }
      result = std::max(result, bound);
    }
    return result + p_rounding;
  };

  // (-T)^j / j!, with absolute error bounds.
  Pair power{1, 0};
  double power_error = 0, weight_sum = 1, weight_lo_sum = 0, weight_error = 0;
  w_hi[0] = 1;
  w_lo[0] = 0;
  for (int j = 1; j <= degree; ++j) {
    const Pair p = two_product(power.hi, -T);
    const double t = power.lo * -T, lo = p.lo + t;
    power_error = T * power_error + u * (std::abs(t) + std::abs(lo));
    power = {p.hi, lo};
    const Pair f = tables.inverse_factorial[j];
    const Pair q = two_product(power.hi, f.hi);
    const double ta = power.hi * f.lo, tb = power.lo * f.hi, tc = ta + tb;
    w_hi[j] = q.hi;
    w_lo[j] = q.lo + tc;
    weight_error +=
        f.hi * power_error +
        u * (std::abs(ta) + std::abs(tb) + std::abs(tc) + std::abs(w_lo[j])) +
        tables.error * std::abs(q.hi) + underflow;
    weight_sum += std::abs(q.hi);
    weight_lo_sum += std::abs(w_lo[j]);
  }

  // Rounding in the accumulation below, per order: |s_hi| <= A and
  // |s_lo| <= AL; each step rounds the two lower-order products, their sum,
  // e, the sum of e and the TwoSum error, and s_lo.
  const double products = weight_sum * max_p_lo + weight_lo_sum * max_p;
  const double exact_lo = u * weight_sum * max_p;
  const double A = 1.01 * (weight_sum + weight_lo_sum) * (max_p + max_p_lo);
  const double AL = 1.01 * (degree * u * A + products + exact_lo);
  const double accumulation =
      1.02 * u * (6 * products + 3 * exact_lo + degree * (u * A + AL));
  const double fixed_error =
      2 * ex * (accumulation + weight_error * (max_p + max_p_lo)) + tail;
  const double budget = options.precision / 4;
  const double total = weight_sum + weight_lo_sum;
  if (!(fixed_error + 2 * total * recurrence_error(false) <= budget) &&
      !(fixed_error + 2 * total * recurrence_error(true) <= budget))
    return false;

  // S_m = sum_j w_j p_{m+j}; independent m lanes vectorize.
  std::copy(p_hi, p_hi + n + 1, s_hi);
  std::copy(p_lo, p_lo + n + 1, s_lo);
  if (fast_fma) {
    for (int j = 1; j <= degree; ++j) {
      const double wh = w_hi[j], wl = w_lo[j];
      const double *ph = p_hi + j, *pl = p_lo + j;
      for (int m = 0; m <= n; ++m) {
        const double p = wh * ph[m];
        const double e = std::fma(wh, ph[m], -p) + (wh * pl[m] + wl * ph[m]);
        const double sum = s_hi[m] + p, t = sum - s_hi[m];
        s_lo[m] += ((s_hi[m] - (sum - t)) + (p - t)) + e;
        s_hi[m] = sum;
      }
    }
  } else {
    // Split each operand once, outside the m loop.
    for (int k = 0; k <= last; ++k) {
      const Pair parts = split(p_hi[k]);
      p_split_hi[k] = parts.hi;
      p_split_lo[k] = parts.lo;
    }
    for (int j = 1; j <= degree; ++j) {
      const double wh = w_hi[j], wl = w_lo[j];
      const Pair weight = split(wh);
      const double *ph = p_hi + j, *pl = p_lo + j;
      const double *qh = p_split_hi + j, *ql = p_split_lo + j;
      for (int m = 0; m <= n; ++m) {
        const double p = wh * ph[m];
        const double e = product_error(weight, {qh[m], ql[m]}, p) +
                         (wh * pl[m] + wl * ph[m]);
        const double sum = s_hi[m] + p, t = sum - s_hi[m];
        s_lo[m] += ((s_hi[m] - (sum - t)) + (p - t)) + e;
        s_hi[m] = sum;
      }
    }
  }
  boys_limit(T, n, w, g);
  for (int m = 0; m <= n; ++m) g[m] = mul_add(-ex, s_hi[m] + s_lo[m], g[m]);
  return true;
}

void propagate_work(double T, int n, Workspace& w, double* g) {
  g[0] = w.g0;
  if (T == 0) {
    for (int m = 1; m <= n; ++m) g[m] = (1 - w.q[m]) / (2. * m + 1);
    return;
  }
  const double e = w.upward_last > 0 || w.seed_order > w.upward_last + 1
                       ? exponential(T, w)
                       : 0;
  for (int m = 0; m < w.upward_last; ++m)
    g[m + 1] = mul_add(m + .5, g[m], .5 * (w.q[m] - e)) / T;
  double next = w.seed;
  if (w.seed_order == n) g[n] = next;
  for (int m = w.seed_order - 1; m > w.upward_last; --m) {
    next = mul_add(2 * T, next, e - w.q[m]) / (2. * m + 1);
    if (m <= n) g[m] = next;
  }
}

void evaluate_work(double R, double T, int n, const Options& options,
                   Workspace& w, double* g) {
  if (R == 0) {
    evaluate_limit(T, n, CutoffLimit::zero, w, g);
    return;
  }
  if (std::isinf(R)) {
    evaluate_limit(T, n, CutoffLimit::boys, w, g);
    return;
  }
  const auto limit = finite_cutoff_limit(R, T, n, options, w);
  if (limit != CutoffLimit::none) {
    evaluate_limit(T, n, limit, w, g);
    return;
  }
  // Below the join at m ~ T - 1/2 the downward branch needs Q_m, which for
  // T >~ 1 only Miller's boundary solve supplies, at several microseconds.
  // The double-word origin expansion gives every order directly; its degree
  // grows with T, which sets the upper end. Small R keeps the radius series.
  if (n > 0 && T > 0 && T < n + .5 && T <= origin_pair_max_t &&
      !(R * R * ((n + 1.5 + T) / 3) <= .25) &&
      origin_pair_taylor(R, T, n, options, w, g))
    return;
  w.reset_request();
  w.g0 = evaluate_g0(R, T, options, w).value;
  if (n == 0) {
    g[0] = w.g0;
    return;
  }
  w.upward_last = 0;
  w.seed_order = -1;
  if (T > 0)
    while (w.upward_last < n && w.upward_last + .5 < T) ++w.upward_last;
  int last = T == 0 ? n : std::max(0, w.upward_last - 1);
  if (T > 0 && w.upward_last < n) {
    const auto seed = direct_seed(R, T, static_cast<unsigned>(n), options, w);
    if (std::isfinite(seed.error)) {
      w.seed = seed.value;
      w.seed_order = n;
      // The join Q_u is unused. If no downward step remains, only the
      // upward inputs Q_0,...,Q_{u-1} are needed (possibly none).
      last = n > w.upward_last + 1 ? n - 1 : w.upward_last - 1;
    } else {
      w.seed = 0;
      w.seed_order =
          zero_seed(T, n, options.precision / 32, options.max_extra_orders);
      last = w.seed_order - 1;
    }
  }
  const double tolerance = T != 0 && last >= 0
                               ? q_tolerances(T, last, options.precision / 2, w)
                               : infinity;
  if (last >= 0) evaluate_q(R, T, last, tolerance, options, w);
  propagate_work(T, n, w, g);
}
}  // namespace

void evaluate(double* Gm, double R, double T, int mmax, double precision,
              Scratch& scratch) {
  const Options options{precision};
  validate(R, T, mmax, options);
  if (!scratch.workspace_)  // moved-from
    scratch.workspace_.reset(new Workspace(mmax));
  evaluate_work(R, T, mmax, options, *scratch.workspace_, Gm);
}
}  // namespace truncated_coulomb
}  // namespace detail
}  // namespace libint2
