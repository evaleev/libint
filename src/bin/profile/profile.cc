/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint compiler.
 *
 *  Libint compiler is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint compiler is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint compiler.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "libint2/util/generated/libint2_params.h"
#ifndef LIBINT2_REALTYPE
#define LIBINT2_REALTYPE double
#endif
#include <cxxabi.h>

#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <numeric>

#include "libint2/boys.h"
#include "libint2/util/timer.h"

#ifndef SKIP_AXPY
#include <mkl_cblas.h>
typedef MKL_INT BLAS_INT;
#endif

#ifdef SKIP_ALL
#define SKIP_AXPY
#define SKIP_DOT
#define SKIP_GEMM
#define SKIP_EXP
#define SKIP_SQRT
#define SKIP_ERF
#define SKIP_CHEBYSHEV
#define SKIP_TAYLOR
#define SKIP_STGNG
#define SKIP_YUKAWA
#endif

#define AVOID_AUTO_VECTORIZATION \
  1  // set to 0 to measure performance within vectorizable loops
     // the default is to measure performance of a single call

using namespace std;

libint2::Timers<1> timer;

enum OperType { f12, f12_2, f12_o_r12, f12_t_f12 };
std::string to_string(OperType o) {
  std::string result;
  switch (o) {
    case f12:
      result = "f12";
      break;
    case f12_2:
      result = "f12_2";
      break;
    case f12_o_r12:
      result = "f12_o_r12";
      break;
    case f12_t_f12:
      result = "f12_t_f12";
      break;
  }
  return result;
}

/* These state variables must be initialized so that they are not all zero. */
namespace dice {
uint32_t x, y, z, w;
void init(unsigned int seed = 159265359) {
  srand(seed);
  x = rand();
  y = rand();
  z = rand();
  w = rand();
}
// uses xorshift128 algorithm, see http://en.wikipedia.org/wiki/Xorshift
uint32_t random_unit32(void) {
  uint32_t t = x ^ (x << 11);
  x = y;
  y = z;
  z = w;
  return w = w ^ (w >> 19) ^ t ^ (t >> 8);
}
};  // namespace dice

/** profiles Kernel by running it repeatedly, and reporting the sum of values
 *
 * @param k the Kernel, needs to have 2 member functions: label and and eval,
 * neither of which takes a parameter
 * @param nrepeats the number of times the kernel is run
 */
template <typename Kernel>
void profile(const Kernel& k, int nrepeats);

template <unsigned InterpolationOrder = 7>
void do_chebyshev(int mmax, int nrepeats);
template <unsigned InterpolationOrder = 7>
void do_taylor(int mmax, int nrepeats);
template <OperType O>
void do_stg6g(int mmax, double T, double rho, int nrepeats);
template <bool exp = false>
void do_stg(int mmax, double T, double U, int nrepeats);

template <typename Real, Real(Function)(Real)>
struct BasicKernel {
  typedef Real ResultType;
  BasicKernel(Real T, std::string label, double T_dec = 0.00001)
      : label_(label), T_(T), T_dec_(T_dec), sum_(0) {}
  std::string label() const { return label_; }
  void eval() const {
    sum_ += Function(T_);
    T_ += T_dec_;
  }
  Real sum() const { return sum_; }
  size_t ops_per_eval() const { return 1; }
  std::string label_;
  mutable Real T_;
  Real T_dec_;
  mutable Real sum_;
};

template <typename Real>
struct VectorOpKernel {
  VectorOpKernel(size_t veclen, size_t nargs, Real T, std::string label)
      : result_(veclen, Real(0)), args_(nargs), label_(label) {
    for (size_t i = 0; i < args_.size(); ++i) {
      args_[i] = new Real[veclen];
      std::fill(args_[i], args_[i] + veclen, T);
    }
  }
  VectorOpKernel(std::initializer_list<size_t> sizes, Real T, std::string label)
      : result_(*sizes.begin(), Real(0)),
        args_(sizes.size() - 1),
        label_(label) {
    auto size_iter = sizes.begin() + 1;
    for (size_t i = 0; i < args_.size(); ++i, ++size_iter) {
      args_[i] = new Real[*size_iter];
      std::fill(args_[i], args_[i] + *size_iter, T);
    }
  }
  ~VectorOpKernel() {
    for (size_t i = 0; i < args_.size(); ++i) {
      delete[] args_[i];
    }
  }
  std::string label() const { return label_; }

  std::vector<Real> result_;
  std::vector<Real*> args_;
  std::string label_;
};

template <typename Real>
struct AXPYKernel : public VectorOpKernel<Real> {
  typedef Real ResultType;
  AXPYKernel(size_t veclen, Real a, Real T, std::string label)
      : VectorOpKernel<Real>(veclen, 1, T, label), a_(a) {
    y_ = &VectorOpKernel<Real>::result_[0];
    x_ = VectorOpKernel<Real>::args_[0];
    n_ = VectorOpKernel<Real>::result_.size();
  }
  void eval() const {
#pragma ivdep
    for (size_t v = 0; v < n_; ++v) y_[v] += a_ * x_[v];
  }
  Real sum() const { return std::accumulate(y_, y_ + n_, Real(0)); }
  size_t ops_per_eval() const { return n_ * 2; }

  Real a_;
  Real* y_;
  Real const* x_;
  BLAS_INT n_;
};

struct DAXPYKernel : public VectorOpKernel<double> {
  typedef double ResultType;
  typedef double Real;
  DAXPYKernel(size_t veclen, Real a, Real T, std::string label)
      : VectorOpKernel<Real>(veclen, 1, T, label), a_(a) {
    y_ = &VectorOpKernel<Real>::result_[0];
    x_ = VectorOpKernel<Real>::args_[0];
    n_ = VectorOpKernel<Real>::result_.size();
  }
  void eval() const { cblas_daxpy(n_, a_, x_, 1, y_, 1); }
  Real sum() const { return std::accumulate(y_, y_ + n_, Real(0)); }
  size_t ops_per_eval() const { return n_ * 2; }

  Real a_;
  Real* y_;
  Real* x_;
  size_t n_;
};

template <typename Real>
struct DOTKernel : public VectorOpKernel<Real> {
  typedef Real ResultType;
  DOTKernel(size_t veclen, Real T, std::string label)
      : VectorOpKernel<Real>(veclen, 2, T, label), result_(0) {
    x1_ = VectorOpKernel<Real>::args_[0];
    x2_ = VectorOpKernel<Real>::args_[1];
    n_ = VectorOpKernel<Real>::result_.size();
  }
  void eval() const {
#pragma ivdep
    for (size_t v = 0; v < n_; ++v) result_ += x1_[v] * x2_[v];
  }
  Real sum() const { return result_; }
  size_t ops_per_eval() const { return n_ * 2; }

  mutable Real result_;
  Real const* x1_;
  Real const* x2_;
  BLAS_INT n_;
};

struct DDOTKernel : public VectorOpKernel<double> {
  typedef double ResultType;
  typedef double Real;
  DDOTKernel(size_t veclen, Real T, std::string label)
      : VectorOpKernel<Real>(veclen, 2, T, label), result_(0) {
    x1_ = VectorOpKernel<Real>::args_[0];
    x2_ = VectorOpKernel<Real>::args_[1];
    n_ = VectorOpKernel<Real>::result_.size();
  }
  void eval() const { result_ += cblas_ddot(n_, x1_, 1, x2_, 1); }
  Real sum() const { return result_; }
  size_t ops_per_eval() const { return n_ * 2; }

  mutable Real result_;
  const Real* x1_;
  const Real* x2_;
  size_t n_;
};

struct DGEMMKernel : public VectorOpKernel<double> {
  typedef double ResultType;
  typedef double Real;
  DGEMMKernel(size_t n, Real T, std::string label)
      : VectorOpKernel<Real>(n * n, 2, T, label), m_(n), n_(n), k_(n) {
    c_ = &VectorOpKernel<Real>::result_[0];
    a_ = VectorOpKernel<Real>::args_[0];
    b_ = VectorOpKernel<Real>::args_[1];
  }
  DGEMMKernel(size_t m, size_t n, size_t k, Real T, std::string label)
      : VectorOpKernel<Real>({m * n, m * k, k * n}, T, label),
        m_(m),
        n_(n),
        k_(k) {
    c_ = &VectorOpKernel<Real>::result_[0];
    a_ = VectorOpKernel<Real>::args_[0];
    b_ = VectorOpKernel<Real>::args_[1];
  }
  void eval() const {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m_, n_, k_, 1.0, a_,
                k_, b_, n_, 1.0, c_, n_);
  }
  Real sum() const { return std::accumulate(c_, c_ + m_ * n_, 0.0); }
  // approximate
  size_t ops_per_eval() const { return m_ * n_ * k_ * 2; }

  Real* c_;        // m by n
  const Real* a_;  // m by k
  const Real* b_;  // k by n
  size_t m_;
  size_t n_;
  size_t k_;
};

const double stg_zeta = 1.0;

int main(int argc, char* argv[]) {
  if (argc < 2 or argc > 4) {
    std::cout << "Description: profiles secondary compute-intensive kernels "
                 "(Boys, daxpy, etc.)"
              << std::endl;
    std::cout << "Usage: profile mmax T nrepeats" << std::endl;
    return 1;
  }
  const int mmax = atoi(argv[1]);
  const double T = atol(argv[2]);
  const double rho = 1.0;
  const int nrepeats = atoi(argv[3]);

  using libint2::simd::VectorSSEDouble;
  const VectorSSEDouble T_sse(T);
#if defined(__AVX__)
  using libint2::simd::VectorAVXDouble;
  const VectorAVXDouble T_avx(T);
#endif

  cout << "mmax = " << mmax << endl;
  cout << " T   = " << T << endl;

  // changes precision of cout to 15 for printing doubles and such.
  cout.precision(15);
  // set the overhead of std::high_resolution_clock; measure carefully
  timer.set_now_overhead(25);

#ifndef SKIP_AXPY
  const int n = 128;
  profile(DAXPYKernel(n, 1.0, 1.0, "daxpy"), nrepeats);
  profile(AXPYKernel<double>(n, 1.0, 1.0, "axpy [double]"), nrepeats);
  profile(AXPYKernel<VectorSSEDouble>(n, 1.0, 1.0, "axpy [SSE]"), nrepeats);
#if defined(__AVX__)
  profile(AXPYKernel<VectorAVXDouble>(n, 1.0, 1.0, "axpy [AVX]"), nrepeats);
#endif
#endif

#ifndef SKIP_DOT
  // const int n = 4096;
  profile(DDOTKernel(n, 1.0, "ddot"), nrepeats);
  profile(DOTKernel<double>(n, 1.0, "dot [double]"), nrepeats);
  profile(DOTKernel<VectorSSEDouble>(n, 1.0, "dot [SSE]"), nrepeats);
#if defined(__AVX__)
  profile(DOTKernel<VectorAVXDouble>(n, 1.0, "dot [AVX]"), nrepeats);
#endif
#endif

#ifndef SKIP_GEMM
  const int nn = 2048;
  profile(DGEMMKernel(nn, 1.0, "dgemm (2048,2048) x (2048,2048)"), 1);
#endif

#ifndef SKIP_EXP
  profile(BasicKernel<double, std::exp>(-T, "exp(-T) [double]", -0.00001),
          nrepeats);
  profile(BasicKernel<VectorSSEDouble, libint2::simd::exp>(-T, "exp(-T) [SSE]",
                                                           -0.00001),
          nrepeats);
#if defined(__AVX__)
  profile(BasicKernel<VectorAVXDouble, libint2::simd::exp>(-T, "exp(-T) [AVX]",
                                                           -0.00001),
          nrepeats);
#endif
#endif
#ifndef SKIP_SQRT
  profile(BasicKernel<double, std::sqrt>(T, "sqrt(T) [double]"), nrepeats);
  profile(BasicKernel<VectorSSEDouble, libint2::simd::sqrt>(T, "sqrt(T) [SSE]"),
          nrepeats);
#if defined(__AVX__)
  profile(BasicKernel<VectorAVXDouble, libint2::simd::sqrt>(T, "sqrt(T) [AVX]"),
          nrepeats);
#endif
#endif
#ifndef SKIP_ERF
  profile(BasicKernel<double, erf>(T, "erf(T) [double]"), nrepeats);
  profile(BasicKernel<VectorSSEDouble, libint2::simd::erf>(T, "erf(T) [SSE]"),
          nrepeats);
#if defined(__AVX__)
  profile(BasicKernel<VectorAVXDouble, libint2::simd::erf>(T, "erf(T) [AVX]"),
          nrepeats);
#endif
  profile(BasicKernel<double, erfc>(T, "erfc(T) [double]"), nrepeats);
  profile(BasicKernel<VectorSSEDouble, libint2::simd::erfc>(T, "erfc(T) [SSE]"),
          nrepeats);
#if defined(__AVX__)
  profile(BasicKernel<VectorAVXDouble, libint2::simd::erfc>(T, "erfc(T) [AVX]"),
          nrepeats);
#endif
#endif
#ifndef SKIP_CHEBYSHEV
  do_chebyshev<7>(mmax, nrepeats);
#endif
#ifndef SKIP_TAYLOR
  //  do_taylor<3>(mmax, nrepeats);
  do_taylor<7>(mmax, nrepeats);
#endif
#ifndef SKIP_STGNG
  do_stg6g<f12>(mmax, T, rho, nrepeats);
  do_stg6g<f12_o_r12>(mmax, T, rho, nrepeats);
  do_stg6g<f12_2>(mmax, T, rho, nrepeats);
  do_stg6g<f12_t_f12>(mmax, T, rho, nrepeats);
#endif
#ifndef SKIP_YUKAWA
  do_stg<true>(mmax, T, rho, nrepeats);
  do_stg<false>(mmax, T, rho, nrepeats);
#endif
  return 0;
}

template <typename Kernel>
void profile(const Kernel& k, int nrepeats) {
  std::cout << "===================== " << k.label()
            << " ======================" << std::endl;

  typedef typename Kernel::ResultType Real;

  timer.clear();
  timer.start(0);

#if AVOID_AUTO_VECTORIZATION
#pragma novector
#endif
  for (int i = 0; i < nrepeats; ++i) {
    k.eval();
  }
  timer.stop(0);

  std::cout << "sum of " << k.label() << ": " << k.sum() << std::endl;

  cout << "Time = " << fixed << timer.read(0) << endl;
  cout << "Rate = " << fixed << nrepeats * k.ops_per_eval() / timer.read(0)
       << endl;
}

template <unsigned InterpolationOrder>
void do_chebyshev(int mmax, int nrepeats) {
  std::cout << "===================== Fm Cheb" << InterpolationOrder
            << " ======================" << std::endl;
  double* Fm_array = new double[mmax + 1];
  double* Fm_array_sum = new double[mmax + 1];
  // std::cout << "alignment of Fm = " << reinterpret_cast<unsigned long
  // int>((char*) Fm_array) % 32ul << " bytes\n";

  // initialize dice
  dice::init();
  const double T_max = 30.0;  // values >= T_max are handled by recursion
  const double scale_unit32_to_T = T_max / std::numeric_limits<uint32_t>::max();

  if (InterpolationOrder != 7) abort();
  typedef libint2::FmEval_Chebyshev7<> fmeval_t;
  fmeval_t fmeval_cheb(mmax);
  std::cout << "done initialization:" << std::endl;
  std::fill(Fm_array_sum, Fm_array_sum + mmax + 1, 0.0);
  timer.clear();
  timer.start(0);
#if AVOID_AUTO_VECTORIZATION
#pragma novector
#endif
  for (int i = 0; i < nrepeats; ++i) {
    const double T = dice::random_unit32() * scale_unit32_to_T;
    // this computes all Fm for up to mmax
    fmeval_cheb.eval(Fm_array, T, mmax);
    for (int m = 0; m <= mmax; ++m) Fm_array_sum[m] += Fm_array[m];
  }
  timer.stop(0);

  std::cout << "sum of Fm:" << std::endl;
  std::copy(Fm_array_sum, Fm_array_sum + mmax + 1,
            std::ostream_iterator<double>(std::cout, "\n"));

  cout << "Time = " << fixed << timer.read(0) << endl;
  cout << "Rate = " << fixed << nrepeats / timer.read(0) << endl;

  delete[] Fm_array;
  delete[] Fm_array_sum;
}

template <unsigned InterpolationOrder>
void do_taylor(int mmax, int nrepeats) {
  std::cout << "===================== Fm Taylor" << InterpolationOrder
            << " ======================" << std::endl;
  double* Fm_array = new double[mmax + 1];
  double* Fm_array_sum = new double[mmax + 1];

  // initialize dice
  dice::init();
  const double T_max = 30.0;  // values >= T_max are handled by recursion
  const double scale_unit32_to_T = T_max / std::numeric_limits<uint32_t>::max();

  libint2::FmEval_Taylor<double, InterpolationOrder> fmeval(mmax);
  std::fill(Fm_array_sum, Fm_array_sum + mmax + 1, 0.0);
  timer.clear();
  timer.start(0);
#if AVOID_AUTO_VECTORIZATION
#pragma novector
#endif
  for (int i = 0; i < nrepeats; ++i) {
    const double T = dice::random_unit32() * scale_unit32_to_T;
    // this computes all Fm for up to mmax
    fmeval.eval(Fm_array, T, mmax);
    for (int m = 0; m <= mmax; ++m) Fm_array_sum[m] += Fm_array[m];
  }
  timer.stop(0);
  std::cout << "sum of Fm (Taylor):" << std::endl;
  std::copy(Fm_array_sum, Fm_array_sum + mmax + 1,
            std::ostream_iterator<double>(std::cout, "\n"));

  cout << "Time = " << fixed << timer.read(0) << endl;
  cout << "Rate = " << fixed << nrepeats / timer.read(0) << endl;

  delete[] Fm_array;
  delete[] Fm_array_sum;
}

template <OperType O>
void do_stg6g(int mmax, double T, double rho, int nrepeats) {
  std::cout << "===================== Gm STG-6G ======================"
            << std::endl;
  double* Gm_array = new double[mmax + 1];
  double* Gm_array_sum = new double[mmax + 1];

  const size_t ng = 6;
  std::vector<std::pair<double, double> > stg_ng(ng);
  // stg_ng[0] = make_pair(4.0001, 1.0);
#if 1
#if HAVE_LAPACK
  libint2::stg_ng_fit(ng, stg_zeta, stg_ng);
#else
  if (stg_zeta != 1.0) {
    throw std::runtime_error("without lapack stg_zeta is hardwired to 0.1");
  }
  stg_ng[0] = make_pair(0.16015391600067220691727771683865433704907890673261,
                        0.20306992259915090794062652264516576964313257462623);
  stg_ng[1] = make_pair(0.58691138376032812074703122125162923674902800850316,
                        0.29474840080158909154305767626309566520662550324323);
  stg_ng[2] = make_pair(1.9052880179050650706871766123674791603725239722860,
                        0.20652315861651088693388092033220845569017370830588);
  stg_ng[3] = make_pair(6.1508186412033182882412135545092215700186355770734,
                        0.13232619560602867340217449542493153700747744735317);
  stg_ng[4] = make_pair(22.558816746266648614747394893787336699960416307706,
                        0.084097701098685716800769730376212853566993914234229);
  stg_ng[5] = make_pair(167.12355778570626548864380047361110482628234458031,
                        0.079234606133959413896805606690618531996594605785539);
#endif
#endif

  std::vector<std::pair<double, double> > stg_ng_sq(stg_ng.size() *
                                                    stg_ng.size());
  for (int b = 0, bk = 0; b < stg_ng.size(); ++b) {
    for (int k = 0; k < stg_ng.size(); ++k, ++bk) {
      const double exp = stg_ng[b].first + stg_ng[k].first;
      const double coef =
          stg_ng[b].second * stg_ng[k].second *
          (O == f12_t_f12 ? 4.0 * stg_ng[b].first * stg_ng[k].first : 1.0);
      stg_ng_sq[bk] = make_pair(exp, coef);
    }
  }

  libint2::GaussianGmEval<
      double, (O == f12 || O == f12_2) ? 0 : (O == f12_o_r12 ? -1 : 2)>
      gtg_eval(mmax, 1e-15);
  std::fill(Gm_array_sum, Gm_array_sum + mmax + 1, 0.0);
  timer.clear();
  timer.start(0);
#if AVOID_AUTO_VECTORIZATION
#pragma novector
#endif
  for (int i = 0; i < nrepeats; ++i) {
    // this computes all Gm for up to mmax
    gtg_eval.eval(Gm_array, rho, T, mmax,
                  (O == f12 || O == f12_o_r12) ? stg_ng : stg_ng_sq);
    for (int m = 0; m <= mmax; ++m) Gm_array_sum[m] += Gm_array[m];
    T += 0.00001;  // to ward-off unrealistic compiler optimizations
    rho += 0.000001;
  }
  timer.stop(0);

  std::cout << "sum of Gm (STG-" << ng << "G O=" << to_string(O)
            << " ):" << std::endl;
  std::copy(Gm_array_sum, Gm_array_sum + mmax + 1,
            std::ostream_iterator<double>(std::cout, "\n"));

  cout << "Time = " << fixed << timer.read(0) << endl;
  cout << "Rate = " << fixed << nrepeats / timer.read(0) << endl;

  delete[] Gm_array;
  delete[] Gm_array_sum;
}

template <bool exp>
void do_stg(int mmax, double T, double rho, int nrepeats) {
  const std::string label = (exp ? "STG" : "Yukawa");
  std::cout << "===================== Gm " << label
            << " ======================" << std::endl;
  double* Gm_array = new double[mmax + 1];
  double* Gm_array_sum = new double[mmax + 1];
  const auto one_over_rho = 1. / rho;

  libint2::TennoGmEval<double> tenno_eval(mmax, 1e-15);
  std::fill(Gm_array_sum, Gm_array_sum + mmax + 1, 0.0);
  timer.clear();
  timer.start(0);

#if AVOID_AUTO_VECTORIZATION
#pragma novector
#endif
  for (int i = 0; i < nrepeats; ++i) {
    // this computes all Gm for up to mmax
    //    asm("#tag1");
    if (exp)
      tenno_eval.eval_slater(Gm_array, one_over_rho, T, mmax, stg_zeta);
    else
      tenno_eval.eval_yukawa(Gm_array, one_over_rho, T, mmax, stg_zeta);
    for (int m = 0; m <= mmax; ++m) Gm_array_sum[m] += Gm_array[m];

    T += 0.00001;  // to ward-off unrealistic compiler optimizations
  }
  timer.stop(0);

  std::cout << "sum of Gm (" << label << "):" << std::endl;
  std::copy(Gm_array_sum, Gm_array_sum + mmax + 1,
            std::ostream_iterator<double>(std::cout, "\n"));

  cout << "Time = " << fixed << timer.read(0) << endl;
  cout << "Rate = " << fixed << nrepeats / timer.read(0) << endl;

  delete[] Gm_array;
  delete[] Gm_array_sum;
}
