/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
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

/// This program tests Libint library by computing 2-body repulsion integrals
/// (4, 3, and 2-center varieties) and (optionally) their derivatives using
/// Libint and a dumb but fool-proof reference method

#include <cmath>
#include <iostream>
#ifdef _MSC_VER
// thanks, https://stackoverflow.com/a/26085827
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <stdint.h>  // portable: uint64_t   MSVC: __int64

// MSVC defines this in winsock2.h!?
typedef struct timeval {
  long tv_sec;
  long tv_usec;
} timeval;

int gettimeofday(struct timeval* tp, struct timezone* tzp) {
  // Note: some broken versions only have 8 trailing zero's, the correct epoch
  // has 9 trailing zero's This magic number is the number of 100 nanosecond
  // intervals since January 1, 1601 (UTC) until 00:00:00 January 1, 1970
  static const uint64_t EPOCH = ((uint64_t)116444736000000000ULL);

  SYSTEMTIME system_time;
  FILETIME file_time;
  uint64_t time;

  GetSystemTime(&system_time);
  SystemTimeToFileTime(&system_time, &file_time);
  time = ((uint64_t)file_time.dwLowDateTime);
  time += ((uint64_t)file_time.dwHighDateTime) << 32;

  tp->tv_sec = (long)((time - EPOCH) / 10000000L);
  tp->tv_usec = (long)(system_time.wMilliseconds * 1000);
  return 0;
}
#else  // _MSC_VER
#include <sys/time.h>
#endif  // _MSC_VER
#include <eri.h>
#include <libint2.h>
#include <libint2/cgshell_ordering.h>
#include <libint2/deriv_iter.h>
#include <libint2/util/memory.h>
#include <prep_libint2.h>

#include <cassert>
#if !LIBINT2_CONSTEXPR_STATICS
#include <libint2/statics_definition.h>
#endif

using namespace std;
using namespace libint2;

const double ABSOLUTE_DEVIATION_THRESHOLD =
    5.0E-14;  // indicate failure if any integral differs in absolute sense by
              // more than this loss of precision in HRR likely limits precision
              // for high-L (e.g. (dp|dd), (dd|dd), etc.)
const double RELATIVE_DEVIATION_THRESHOLD =
    1.0E-9;  // indicate failure if any integral differs in relative sense by
             // more than this

double absolute_derivation_threshold(int deriv_order) {
  return ABSOLUTE_DEVIATION_THRESHOLD *
         std::pow(3., deriv_order > 2 ? deriv_order - 2 : 0);
}

/// change to true to skip verification and do some timing simulation
const bool do_timing_only = false;

libint2::FmEval_Chebyshev7<double> fmeval_chebyshev(std::max(LIBINT_MAX_AM, 4) *
                                                        4 +
                                                    2);
libint2::FmEval_Taylor<double, 6> fmeval_taylor(std::max(LIBINT_MAX_AM, 4) * 4 +
                                                    2,
                                                1e-15);

namespace {

std::string am_to_symbol(unsigned int l, bool contracted = false) {
  std::string result;
  const size_t lmax_plus_1 = sizeof(LIBINT_AM2SYMBOL) - 1;
  do {
    const unsigned int digit = l % lmax_plus_1;
    char letter = LIBINT_AM2SYMBOL[digit];
    if (contracted) letter = toupper(letter);
    result.insert(result.begin(), letter);
    l /= lmax_plus_1;
  } while (l != 0);

  return result;
}
}  // namespace

// test 4, 3, and 2-center integrals
#ifdef INCLUDE_ERI
bool test_4eri(unsigned int deriv_order, unsigned int lmax_max);
#endif
#ifdef INCLUDE_ERI3
bool test_3eri(unsigned int deriv_order, unsigned int lmax_max);
#endif
#ifdef INCLUDE_ERI2
bool test_2eri(unsigned int deriv_order, unsigned int lmax_max);
#endif

/// give optional derivative order (default = 0, i.e. regular integrals)
int main(int argc, char** argv) {
  assert(argc == 1 || argc == 2 || argc == 3);
  const unsigned int deriv_order =
      (argc == 2 || argc == 3) ? atoi(argv[1]) : 0u;
  const unsigned int lmax_max = (argc == 3) ? atoi(argv[2]) : UINT_MAX;

  // static initialization of the library (only needs to happen once per
  // process)
  LIBINT2_PREFIXED_NAME(libint2_static_init)();

  bool success = true;
  // run the tests
#ifdef INCLUDE_ERI
  if (deriv_order <= 4) success &= test_4eri(deriv_order, lmax_max);
#endif
#ifdef INCLUDE_ERI3
  if (deriv_order <= 4) success &= test_3eri(deriv_order, lmax_max);
#endif
#ifdef INCLUDE_ERI2
  if (deriv_order <= 4) success &= test_2eri(deriv_order, lmax_max);
#endif

  // cleanup static library data (once per process)
  LIBINT2_PREFIXED_NAME(libint2_static_cleanup)();

  return success ? 0 : 1;
}

#ifdef INCLUDE_ERI
bool test_4eri(unsigned int deriv_order, unsigned int lmax_max) {
  if (deriv_order > INCLUDE_ERI) return true;

  bool test_success = true;

  // record start wall time
  struct timeval tod;
  gettimeofday(&tod, 0);
  const double start_wall_time = tod.tv_sec + 0.000001 * tod.tv_usec;

  typedef unsigned int uint;

  const uint veclen = LIBINT2_MAX_VECLEN;
  const uint max_contrdepth = 3;
  const uint max_contrdepth4 =
      max_contrdepth * max_contrdepth * max_contrdepth * max_contrdepth;

  unsigned int lmax;
  if (deriv_order == 0) lmax = LIBINT2_MAX_AM_eri;
#if INCLUDE_ERI >= 1
  if (deriv_order == 1) lmax = LIBINT2_MAX_AM_eri1;
#endif
#if INCLUDE_ERI >= 2
  if (deriv_order == 2) lmax = LIBINT2_MAX_AM_eri2;
#endif
#if INCLUDE_ERI >= 3
  if (deriv_order == 3) lmax = LIBINT2_MAX_AM_eri3;
#endif
#if INCLUDE_ERI >= 4
  if (deriv_order == 4) lmax = LIBINT2_MAX_AM_eri4;
#endif
  Libint_t* inteval = libint2::malloc<Libint_t>(max_contrdepth4);
  if (deriv_order == 0)
    LIBINT2_PREFIXED_NAME(libint2_init_eri)(&inteval[0], lmax, 0);
#if INCLUDE_ERI >= 1
  if (deriv_order == 1)
    LIBINT2_PREFIXED_NAME(libint2_init_eri1)(&inteval[0], lmax, 0);
#endif
#if INCLUDE_ERI >= 2
  if (deriv_order == 2)
    LIBINT2_PREFIXED_NAME(libint2_init_eri2)(&inteval[0], lmax, 0);
#endif
#if INCLUDE_ERI >= 3
  if (deriv_order == 3)
    LIBINT2_PREFIXED_NAME(libint2_init_eri3)(&inteval[0], lmax, 0);
#endif
#if INCLUDE_ERI >= 4
  if (deriv_order == 4)
    LIBINT2_PREFIXED_NAME(libint2_init_eri4)(&inteval[0], lmax, 0);
#endif
#ifdef LIBINT2_FLOP_COUNT
  LIBINT2_PREFIXED_NAME(libint2_init_flopcounter)(&inteval[0], max_contrdepth4);
#endif

  lmax = std::min(lmax_max, lmax);

  for (unsigned int l0 = 0; l0 <= lmax; ++l0) {
    for (unsigned int l1 = 0; l1 <= lmax; ++l1) {
      for (unsigned int l2 = 0; l2 <= lmax; ++l2) {
        for (unsigned int l3 = 0; l3 <= lmax; ++l3) {
          // record start wall time
          struct timeval tod;
          gettimeofday(&tod, 0);
          const double start_wall_time = tod.tv_sec + 0.000001 * tod.tv_usec;

#if LIBINT_CONTRACTED_INTS
          const uint contrdepth =
              do_timing_only
                  ? std::min((4 * lmax + 4) / (l0 + l1 + l2 + l3 + 4),
                             max_contrdepth)
                  : max_contrdepth;
#else
          const uint contrdepth = 1;
#endif
          const uint contrdepth4 =
              contrdepth * contrdepth * contrdepth * contrdepth;

          // can compute this? skip, if not.
          // there are many reason why Libint could not compute a type of
          // integrals for example, Libint does not compute (ss|ss) integrals
          // (although it does compute derivatives of (ss|ss) another reason is
          // a given integral type is not unique and can be computed using other
          // functions in Libint
          if (deriv_order == 0 &&
              LIBINT2_PREFIXED_NAME(libint2_build_eri)[l0][l1][l2][l3] == 0)
            continue;
#if INCLUDE_ERI >= 1
          if (deriv_order == 1 &&
              LIBINT2_PREFIXED_NAME(libint2_build_eri1)[l0][l1][l2][l3] == 0)
            continue;
#endif
#if INCLUDE_ERI >= 2
          if (deriv_order == 2 &&
              LIBINT2_PREFIXED_NAME(libint2_build_eri2)[l0][l1][l2][l3] == 0)
            continue;
#endif
#if INCLUDE_ERI >= 3
          if (deriv_order == 3 &&
              LIBINT2_PREFIXED_NAME(libint2_build_eri3)[l0][l1][l2][l3] == 0)
            continue;
#endif
#if INCLUDE_ERI >= 4
          if (deriv_order == 4 &&
              LIBINT2_PREFIXED_NAME(libint2_build_eri4)[l0][l1][l2][l3] == 0)
            continue;
#endif

          unsigned int am[4];
          am[0] = l0;
          am[1] = l1;
          am[2] = l2;
          am[3] = l3;
          RandomShellSet<4u> rsqset(am, veclen, contrdepth, 172530558);

          CartesianDerivIterator<4> diter(deriv_order);
          const unsigned int nderiv = diter.range_size();

          const double* A = &(rsqset.R[0][0]);
          const double* B = &(rsqset.R[1][0]);
          const double* C = &(rsqset.R[2][0]);
          const double* D = &(rsqset.R[3][0]);
          LIBINT2_REF_REALTYPE Aref[3];
          for (int i = 0; i < 3; ++i) Aref[i] = A[i];
          LIBINT2_REF_REALTYPE Bref[3];
          for (int i = 0; i < 3; ++i) Bref[i] = B[i];
          LIBINT2_REF_REALTYPE Cref[3];
          for (int i = 0; i < 3; ++i) Cref[i] = C[i];
          LIBINT2_REF_REALTYPE Dref[3];
          for (int i = 0; i < 3; ++i) Dref[i] = D[i];

          const int nrepeats = do_timing_only
                                   ? 50 * (lmax - l0 + 1) * (lmax - l1 + 1) *
                                         (lmax - l2 + 1) * (lmax - l3 + 1)
                                   : 1;

          cout << (do_timing_only ? "Timing " : "Testing ") << " ("
               << am_to_symbol(am[0]) << am_to_symbol(am[1]) << "|"
               << am_to_symbol(am[2]) << am_to_symbol(am[3]) << ") ";
          if (deriv_order > 0) {
            cout << " deriv order = " << deriv_order;
          }
          if (do_timing_only) {
            cout << " contrdepth = " << contrdepth
                 << " #(repeats) = " << nrepeats;
          }
          cout << ": ";

          for (int k = 0; k < nrepeats; ++k) {
            // this prepares the data
            prep_libint2(&inteval[0], rsqset, 0, deriv_order);

            //  now use Libint to compute
            double scale_target = 1.0;
#if LIBINT_ACCUM_INTS
            // if accumulating integrals, zero out first, then compute twice
            inteval[0].zero_out_targets = 1;
            scale_target = 0.5;
            if (deriv_order == 0) LIBINT2_PREFIXED_NAME(libint2_build_eri)
            [am[0]][am[1]][am[2]][am[3]](&inteval[0]);
#if INCLUDE_ERI >= 1
            else if (deriv_order == 1) LIBINT2_PREFIXED_NAME(
                libint2_build_eri1)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);
#endif
#if INCLUDE_ERI >= 2
            else if (deriv_order == 2) LIBINT2_PREFIXED_NAME(
                libint2_build_eri2)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);
#endif
#if INCLUDE_ERI >= 3
            else if (deriv_order == 3) LIBINT2_PREFIXED_NAME(
                libint2_build_eri3)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);
#endif
#if INCLUDE_ERI >= 4
            else if (deriv_order == 4) LIBINT2_PREFIXED_NAME(
                libint2_build_eri4)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);
#endif
#endif
#if LIBINT_CONTRACTED_INTS
            inteval[0].contrdepth = contrdepth4;
#endif
            if (deriv_order == 0) LIBINT2_PREFIXED_NAME(libint2_build_eri)
            [am[0]][am[1]][am[2]][am[3]](&inteval[0]);
#if INCLUDE_ERI >= 1
            else if (deriv_order == 1) LIBINT2_PREFIXED_NAME(
                libint2_build_eri1)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);
#endif
#if INCLUDE_ERI >= 2
            else if (deriv_order == 2) LIBINT2_PREFIXED_NAME(
                libint2_build_eri2)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);
#endif
#if INCLUDE_ERI >= 3
            else if (deriv_order == 3) LIBINT2_PREFIXED_NAME(
                libint2_build_eri3)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);
#endif
#if INCLUDE_ERI >= 4
            else if (deriv_order == 4) LIBINT2_PREFIXED_NAME(
                libint2_build_eri4)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);
#endif

            if (not do_timing_only) {
              // compare Libint integrals against the reference method
              // since the reference implementation computes integrals one at a
              // time (not one shell-set at a time) the outer loop is over the
              // basis functions

              bool success = true;
              int ijkl = 0;

              int l0, m0, n0;
              FOR_CART(l0, m0, n0, am[0])

              int l1, m1, n1;
              FOR_CART(l1, m1, n1, am[1])

              int l2, m2, n2;
              FOR_CART(l2, m2, n2, am[2])

              int l3, m3, n3;
              FOR_CART(l3, m3, n3, am[3])

              for (uint v = 0; v < veclen; v++) {
                //
                // compute reference integrals
                //
                std::vector<LIBINT2_REF_REALTYPE> ref_eri(nderiv, 0.0);

                uint p0123 = 0;
                for (uint p0 = 0; p0 < contrdepth; p0++) {
                  for (uint p1 = 0; p1 < contrdepth; p1++) {
                    for (uint p2 = 0; p2 < contrdepth; p2++) {
                      for (uint p3 = 0; p3 < contrdepth; p3++, p0123++) {
                        const LIBINT2_REF_REALTYPE alpha0 =
                            rsqset.exp[0][v][p0];
                        const LIBINT2_REF_REALTYPE alpha1 =
                            rsqset.exp[1][v][p1];
                        const LIBINT2_REF_REALTYPE alpha2 =
                            rsqset.exp[2][v][p2];
                        const LIBINT2_REF_REALTYPE alpha3 =
                            rsqset.exp[3][v][p3];

                        const LIBINT2_REF_REALTYPE c0 = rsqset.coef[0][v][p0];
                        const LIBINT2_REF_REALTYPE c1 = rsqset.coef[1][v][p1];
                        const LIBINT2_REF_REALTYPE c2 = rsqset.coef[2][v][p2];
                        const LIBINT2_REF_REALTYPE c3 = rsqset.coef[3][v][p3];
                        const LIBINT2_REF_REALTYPE c0123 = c0 * c1 * c2 * c3;

                        CartesianDerivIterator<4> diter(deriv_order);
                        bool last_deriv = false;
                        unsigned int di = 0;
                        do {
                          ref_eri[di++] +=
                              c0123 * eri(&(*diter)[0], l0, m0, n0, alpha0,
                                          Aref, l1, m1, n1, alpha1, Bref, l2,
                                          m2, n2, alpha2, Cref, l3, m3, n3,
                                          alpha3, Dref, 0);
                          last_deriv = diter.last();
                          if (!last_deriv) diter.next();
                        } while (!last_deriv);
                      }
                    }
                  }
                }

                //
                // extract Libint integrals
                //
                std::vector<LIBINT2_REALTYPE> new_eri;
                for (auto d = 0; d != nderiv; ++d)
                  new_eri.push_back(scale_target *
                                    inteval[0].targets[d][ijkl * veclen + v]);

                //
                // compare reference and libint integrals
                //
                for (unsigned int di = 0; di < nderiv; ++di) {
                  const LIBINT2_REF_REALTYPE abs_error =
                      abs(ref_eri[di] - LIBINT2_REF_REALTYPE(new_eri[di]));
                  const LIBINT2_REF_REALTYPE relabs_error =
                      abs(abs_error / ref_eri[di]);
                  if (relabs_error > RELATIVE_DEVIATION_THRESHOLD &&
                      abs_error > absolute_derivation_threshold(deriv_order)) {
                    std::cout << "Elem " << ijkl << " di= " << di << " v=" << v
                              << " : ref = " << ref_eri[di]
                              << " libint = " << new_eri[di]
                              << " relabs_error = " << relabs_error
                              << " abs_error = " << abs_error << endl;
                    success = false;
                  }
                }

              }  // end of vector loop
              ++ijkl;
              END_FOR_CART
              END_FOR_CART
              END_FOR_CART
              END_FOR_CART

              cout << (success ? "ok" : "failed") << std::endl;

              test_success &= success;
            }  // checking computed values vs. the reference

          }  // end of nrepeats

          if (do_timing_only) {
            // record end wall time, compute total wall time spent here
            gettimeofday(&tod, 0);
            const double end_wall_time = tod.tv_sec + 0.000001 * tod.tv_usec;
            std::cout << "wall time = " << (end_wall_time - start_wall_time)
                      << " seconds" << std::endl;
          }
        }
      }
    }
  }

  if (deriv_order == 0) LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);
#if INCLUDE_ERI >= 1
  if (deriv_order == 1)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri1)(&inteval[0]);
#endif
#if INCLUDE_ERI >= 2
  if (deriv_order == 2)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri2)(&inteval[0]);
#endif
#if INCLUDE_ERI >= 3
  if (deriv_order == 3)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri3)(&inteval[0]);
#endif
#if INCLUDE_ERI >= 4
  if (deriv_order == 4)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri4)(&inteval[0]);
#endif
  free(inteval);

  // record end wall time, compute total wall time spent here
  gettimeofday(&tod, 0);
  const double end_wall_time = tod.tv_sec + 0.000001 * tod.tv_usec;
  std::cout << "wall time = " << (end_wall_time - start_wall_time) << " seconds"
            << std::endl;

  return test_success;
}
#endif  // INCLUDE_ERI

#ifdef INCLUDE_ERI3
bool test_3eri(unsigned int deriv_order, unsigned int lmax_max) {
  if (deriv_order > INCLUDE_ERI3) return true;

  bool test_success = true;

  // record start wall time
  struct timeval tod;
  gettimeofday(&tod, 0);
  const double start_wall_time = tod.tv_sec + 0.000001 * tod.tv_usec;

  typedef unsigned int uint;

  const uint veclen = LIBINT2_MAX_VECLEN;
  const uint max_contrdepth = 3;
  const uint max_contrdepth3 = max_contrdepth * max_contrdepth * max_contrdepth;

  unsigned int lmax_default = LIBINT2_MAX_AM;
#if defined(LIBINT2_MAX_AM1)
  if (deriv_order == 1) lmax_default = LIBINT2_MAX_AM1;
#endif
#if defined(LIBINT2_MAX_AM2)
  if (deriv_order == 2) lmax_default = LIBINT2_MAX_AM2;
#endif
#if defined(LIBINT2_MAX_AM3)
  if (deriv_order == 3) lmax_default = LIBINT2_MAX_AM3;
#endif
#if defined(LIBINT2_MAX_AM4)
  if (deriv_order == 4) lmax_default = LIBINT2_MAX_AM4;
#endif

  unsigned int lmax;
  if (deriv_order == 0) lmax = LIBINT2_MAX_AM_3eri;
#if INCLUDE_ERI3 >= 1
  if (deriv_order == 1) lmax = LIBINT2_MAX_AM_3eri1;
#endif
#if INCLUDE_ERI3 >= 2
  if (deriv_order == 2) lmax = LIBINT2_MAX_AM_3eri2;
#endif
#if INCLUDE_ERI3 >= 3
  if (deriv_order == 3) lmax = LIBINT2_MAX_AM_3eri3;
#endif
#if INCLUDE_ERI3 >= 4
  if (deriv_order == 4) lmax = LIBINT2_MAX_AM_3eri4;
#endif

  const auto lmax_init = std::max(lmax, lmax_default);

  Libint_t* inteval = libint2::malloc<Libint_t>(max_contrdepth3);
  if (deriv_order == 0)
    LIBINT2_PREFIXED_NAME(libint2_init_3eri)(&inteval[0], lmax_init, 0);
#if INCLUDE_ERI3 >= 1
  if (deriv_order == 1)
    LIBINT2_PREFIXED_NAME(libint2_init_3eri1)(&inteval[0], lmax_init, 0);
#endif
#if INCLUDE_ERI3 >= 2
  if (deriv_order == 2)
    LIBINT2_PREFIXED_NAME(libint2_init_3eri2)(&inteval[0], lmax_init, 0);
#endif
#if INCLUDE_ERI3 >= 3
  if (deriv_order == 3)
    LIBINT2_PREFIXED_NAME(libint2_init_3eri3)(&inteval[0], lmax_init, 0);
#endif
#if INCLUDE_ERI3 >= 4
  if (deriv_order == 4)
    LIBINT2_PREFIXED_NAME(libint2_init_3eri4)(&inteval[0], lmax_init, 0);
#endif
#ifdef LIBINT2_FLOP_COUNT
  LIBINT2_PREFIXED_NAME(libint2_init_flopcounter)(&inteval[0], max_contrdepth3);
#endif

  auto lmax0 = std::min(lmax_max, lmax);
  auto lmax1 = std::min(lmax_max, lmax_default);
  auto lmax2 = std::min(lmax_max, lmax_default);

  // reference code only supports Cartesian gaussians
#if ERI3_PURE_SH
  lmax0 = std::min(lmax0, 1u);
#endif

  for (unsigned int l0 = 0; l0 <= lmax0; ++l0) {
    for (unsigned int l1 = 0; l1 <= lmax1; ++l1) {
      for (unsigned int l2 = 0; l2 <= lmax2; ++l2) {
        // record start wall time
        struct timeval tod;
        gettimeofday(&tod, 0);
        const double start_wall_time = tod.tv_sec + 0.000001 * tod.tv_usec;

        // can compute this? skip, if not
        if (deriv_order == 0 &&
            LIBINT2_PREFIXED_NAME(libint2_build_3eri)[l0][l1][l2] == 0)
          continue;
#if INCLUDE_ERI3 >= 1
        if (deriv_order == 1 &&
            LIBINT2_PREFIXED_NAME(libint2_build_3eri1)[l0][l1][l2] == 0)
          continue;
#endif
#if INCLUDE_ERI3 >= 2
        if (deriv_order == 2 &&
            LIBINT2_PREFIXED_NAME(libint2_build_3eri2)[l0][l1][l2] == 0)
          continue;
#endif
#if INCLUDE_ERI3 >= 3
        if (deriv_order == 3 &&
            LIBINT2_PREFIXED_NAME(libint2_build_3eri3)[l0][l1][l2] == 0)
          continue;
#endif
#if INCLUDE_ERI3 >= 4
        if (deriv_order == 4 &&
            LIBINT2_PREFIXED_NAME(libint2_build_3eri4)[l0][l1][l2] == 0)
          continue;
#endif

#if LIBINT_CONTRACTED_INTS
        const uint contrdepth =
            do_timing_only
                ? std::min((3 * lmax + 3) / (l0 + l1 + l2 + 3), max_contrdepth)
                : max_contrdepth;
#else
        const uint contrdepth = 1;
#endif
        const uint contrdepth3 = contrdepth * contrdepth * contrdepth;

        unsigned int am[3];
        am[0] = l0;
        am[1] = l1;
        am[2] = l2;
        RandomShellSet<3u> rsqset(am, veclen, contrdepth, 201375807);

        CartesianDerivIterator<3> diter(deriv_order);
        const unsigned int nderiv = diter.range_size();

        const double* A = &(rsqset.R[0][0]);
        const double* B = &(rsqset.R[1][0]);
        const double* C = &(rsqset.R[2][0]);
        LIBINT2_REF_REALTYPE Aref[3];
        for (int i = 0; i < 3; ++i) Aref[i] = A[i];
        LIBINT2_REF_REALTYPE Bref[3];
        for (int i = 0; i < 3; ++i) Bref[i] = B[i];
        LIBINT2_REF_REALTYPE Cref[3];
        for (int i = 0; i < 3; ++i) Cref[i] = C[i];

        const int nrepeats =
            do_timing_only
                ? 1000 * (lmax - l0 + 1) * (lmax - l1 + 1) * (lmax - l2 + 1)
                : 1;

        cout << (do_timing_only ? "Timing " : "Testing ") << "("
             << am_to_symbol(am[0]) << "|" << am_to_symbol(am[1])
             << am_to_symbol(am[2]) << ") ";
        if (deriv_order > 0) {
          cout << " deriv order = " << deriv_order;
        }
        if (do_timing_only) {
          cout << " contrdepth = " << contrdepth
               << " #(repeats) = " << nrepeats;
        }
        cout << ": ";

        for (int k = 0; k < nrepeats; ++k) {
          prep_libint2(inteval, rsqset, 0, deriv_order);

          double scale_target = 1.0;
#if LIBINT_ACCUM_INTS
          // if accumulating integrals, zero out first, then compute twice
          inteval[0].zero_out_targets = 1;
          scale_target = 0.5;
          if (deriv_order == 0) LIBINT2_PREFIXED_NAME(libint2_build_3eri)
          [am[0]][am[1]][am[2]](&inteval[0]);
#if INCLUDE_ERI3 >= 1
          if (deriv_order == 1) LIBINT2_PREFIXED_NAME(libint2_build_3eri1)
          [am[0]][am[1]][am[2]](&inteval[0]);
#endif
#if INCLUDE_ERI3 >= 2
          if (deriv_order == 2) LIBINT2_PREFIXED_NAME(libint2_build_3eri2)
          [am[0]][am[1]][am[2]](&inteval[0]);
#endif
#if INCLUDE_ERI3 >= 3
          if (deriv_order == 3) LIBINT2_PREFIXED_NAME(libint2_build_3eri3)
          [am[0]][am[1]][am[2]](&inteval[0]);
#endif
#if INCLUDE_ERI3 >= 4
          if (deriv_order == 4) LIBINT2_PREFIXED_NAME(libint2_build_3eri4)
          [am[0]][am[1]][am[2]](&inteval[0]);
#endif
#endif
#if LIBINT_CONTRACTED_INTS
          inteval[0].contrdepth = contrdepth3;
#endif
          if (deriv_order == 0) LIBINT2_PREFIXED_NAME(libint2_build_3eri)
          [am[0]][am[1]][am[2]](&inteval[0]);
#if INCLUDE_ERI3 >= 1
          if (deriv_order == 1) LIBINT2_PREFIXED_NAME(libint2_build_3eri1)
          [am[0]][am[1]][am[2]](&inteval[0]);
#endif
#if INCLUDE_ERI3 >= 2
          if (deriv_order == 2) LIBINT2_PREFIXED_NAME(libint2_build_3eri2)
          [am[0]][am[1]][am[2]](&inteval[0]);
#endif
#if INCLUDE_ERI3 >= 3
          if (deriv_order == 3) LIBINT2_PREFIXED_NAME(libint2_build_3eri3)
          [am[0]][am[1]][am[2]](&inteval[0]);
#endif
#if INCLUDE_ERI3 >= 4
          if (deriv_order == 4) LIBINT2_PREFIXED_NAME(libint2_build_3eri4)
          [am[0]][am[1]][am[2]](&inteval[0]);
#endif

          if (not do_timing_only) {
            bool success = true;
            int ijk = 0;

            int l0, m0, n0;
            FOR_CART(l0, m0, n0, am[0])

            int l1, m1, n1;
            FOR_CART(l1, m1, n1, am[1])

            int l2, m2, n2;
            FOR_CART(l2, m2, n2, am[2])

            for (uint v = 0; v < veclen; v++) {
              std::vector<LIBINT2_REF_REALTYPE> ref_eri(nderiv, 0.0);

              uint p012 = 0;
              for (uint p0 = 0; p0 < contrdepth; p0++) {
                for (uint p1 = 0; p1 < contrdepth; p1++) {
                  for (uint p2 = 0; p2 < contrdepth; p2++, p012++) {
                    const LIBINT2_REF_REALTYPE alpha0 = rsqset.exp[0][v][p0];
                    const LIBINT2_REF_REALTYPE alpha1 = rsqset.exp[1][v][p1];
                    const LIBINT2_REF_REALTYPE alpha2 = rsqset.exp[2][v][p2];

                    const LIBINT2_REF_REALTYPE c0 = rsqset.coef[0][v][p0];
                    const LIBINT2_REF_REALTYPE c1 = rsqset.coef[1][v][p1];
                    const LIBINT2_REF_REALTYPE c2 = rsqset.coef[2][v][p2];
                    const LIBINT2_REF_REALTYPE c012 = c0 * c1 * c2;

                    CartesianDerivIterator<3> diter(deriv_order);
                    bool last_deriv = false;
                    unsigned int di = 0;
                    do {
                      // convert 3-center deriv indices into 4-center deriv
                      // indices
                      unsigned int deriv_level[12];
                      auto* deriv3_level = &(*diter)[0];
                      std::copy(deriv3_level, deriv3_level + 3, deriv_level);
                      std::fill(deriv_level + 3, deriv_level + 6, 0u);
                      std::copy(deriv3_level + 3, deriv3_level + 9,
                                deriv_level + 6);

                      ref_eri[di++] +=
                          c012 * eri(deriv_level, l0, m0, n0, alpha0, Aref, 0u,
                                     0u, 0u, 0.0, Aref, l1, m1, n1, alpha1,
                                     Bref, l2, m2, n2, alpha2, Cref, 0);
                      last_deriv = diter.last();
                      if (!last_deriv) diter.next();
                    } while (!last_deriv);
                  }
                }
              }

              //
              // extract Libint integrals
              //
              std::vector<LIBINT2_REALTYPE> new_eri;
              for (auto d = 0; d != nderiv; ++d)
                new_eri.push_back(scale_target *
                                  inteval[0].targets[d][ijk * veclen + v]);

              for (unsigned int di = 0; di < nderiv; ++di) {
                const LIBINT2_REF_REALTYPE abs_error =
                    abs(ref_eri[di] - LIBINT2_REF_REALTYPE(new_eri[di]));
                const LIBINT2_REF_REALTYPE relabs_error =
                    abs(abs_error / ref_eri[di]);
                if (relabs_error > RELATIVE_DEVIATION_THRESHOLD &&
                    abs_error > absolute_derivation_threshold(deriv_order)) {
                  std::cout << "Elem " << ijk << " di= " << di << " v=" << v
                            << " : ref = " << ref_eri[di]
                            << " libint = " << new_eri[di]
                            << " relabs_error = " << relabs_error << endl;
                  success = false;
                }
              }

            }  // end of vector loop
            ++ijk;
            END_FOR_CART
            END_FOR_CART
            END_FOR_CART

            cout << (success ? "ok" : "failed") << endl;

            test_success &= success;
          }  // checking computed values vs. the reference

        }  // end of nrepeats

        if (do_timing_only) {
          // record end wall time, compute total wall time spent here
          gettimeofday(&tod, 0);
          const double end_wall_time = tod.tv_sec + 0.000001 * tod.tv_usec;
          std::cout << "wall time = " << (end_wall_time - start_wall_time)
                    << " seconds" << std::endl;
        }
      }
    }
  }

  if (deriv_order == 0)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_3eri)(&inteval[0]);
#if INCLUDE_ERI3 >= 1
  if (deriv_order == 1)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_3eri1)(&inteval[0]);
#endif
#if INCLUDE_ERI3 >= 2
  if (deriv_order == 2)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_3eri2)(&inteval[0]);
#endif
#if INCLUDE_ERI3 >= 3
  if (deriv_order == 3)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_3eri3)(&inteval[0]);
#endif
#if INCLUDE_ERI3 >= 4
  if (deriv_order == 4)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_3eri4)(&inteval[0]);
#endif
  free(inteval);

  // record end wall time, compute total wall time spent here
  gettimeofday(&tod, 0);
  const double end_wall_time = tod.tv_sec + 0.000001 * tod.tv_usec;
  std::cout << "wall time = " << (end_wall_time - start_wall_time) << " seconds"
            << std::endl;

  return test_success;
}
#endif  // INCLUDE_ERI3

#ifdef INCLUDE_ERI2
bool test_2eri(unsigned int deriv_order, unsigned int lmax_max) {
  if (deriv_order > INCLUDE_ERI2) return true;

  bool test_success = true;

  typedef unsigned int uint;

  const uint veclen = LIBINT2_MAX_VECLEN;
#if LIBINT_CONTRACTED_INTS
  const uint contrdepth = 3;
#else
  const uint contrdepth = 1;
#endif
  const uint contrdepth2 = contrdepth * contrdepth;

  unsigned int lmax;
  if (deriv_order == 0) lmax = LIBINT2_MAX_AM_2eri;
#if INCLUDE_ERI2 >= 1
  if (deriv_order == 1) lmax = LIBINT2_MAX_AM_2eri1;
#endif
#if INCLUDE_ERI2 >= 2
  if (deriv_order == 2) lmax = LIBINT2_MAX_AM_2eri2;
#endif
#if INCLUDE_ERI2 >= 3
  if (deriv_order == 3) lmax = LIBINT2_MAX_AM_2eri3;
#endif
#if INCLUDE_ERI2 >= 4
  if (deriv_order == 4) lmax = LIBINT2_MAX_AM_2eri4;
#endif
  Libint_t* inteval = libint2::malloc<Libint_t>(contrdepth2);
  if (deriv_order == 0)
    LIBINT2_PREFIXED_NAME(libint2_init_2eri)(&inteval[0], lmax, 0);
#if INCLUDE_ERI2 >= 1
  if (deriv_order == 1)
    LIBINT2_PREFIXED_NAME(libint2_init_2eri1)(&inteval[0], lmax, 0);
#endif
#if INCLUDE_ERI2 >= 2
  if (deriv_order == 2)
    LIBINT2_PREFIXED_NAME(libint2_init_2eri2)(&inteval[0], lmax, 0);
#endif
#if INCLUDE_ERI2 >= 3
  if (deriv_order == 3)
    LIBINT2_PREFIXED_NAME(libint2_init_2eri3)(&inteval[0], lmax, 0);
#endif
#if INCLUDE_ERI2 >= 4
  if (deriv_order == 4)
    LIBINT2_PREFIXED_NAME(libint2_init_2eri4)(&inteval[0], lmax, 0);
#endif
#ifdef LIBINT2_FLOP_COUNT
  LIBINT2_PREFIXED_NAME(libint2_init_flopcounter)(&inteval[0], contrdepth2);
#endif

  lmax = std::min(lmax_max, lmax);

  // reference code only supports Cartesian gaussians
#if ERI2_PURE_SH
  lmax = std::min(lmax, 1u);
#endif

  for (unsigned int l0 = 0; l0 <= lmax; ++l0) {
    for (unsigned int l1 = 0; l1 <= lmax; ++l1) {
      // can compute this? skip, if not
      if (deriv_order == 0 &&
          LIBINT2_PREFIXED_NAME(libint2_build_2eri)[l0][l1] == 0)
        continue;
#if INCLUDE_ERI2 >= 1
      if (deriv_order == 1 &&
          LIBINT2_PREFIXED_NAME(libint2_build_2eri1)[l0][l1] == 0)
        continue;
#endif
#if INCLUDE_ERI2 >= 2
      if (deriv_order == 2 &&
          LIBINT2_PREFIXED_NAME(libint2_build_2eri2)[l0][l1] == 0)
        continue;
#endif
#if INCLUDE_ERI2 >= 3
      if (deriv_order == 3 &&
          LIBINT2_PREFIXED_NAME(libint2_build_2eri3)[l0][l1] == 0)
        continue;
#endif
#if INCLUDE_ERI2 >= 4
      if (deriv_order == 4 &&
          LIBINT2_PREFIXED_NAME(libint2_build_2eri4)[l0][l1] == 0)
        continue;
#endif

      unsigned int am[2];
      am[0] = l0;
      am[1] = l1;
      RandomShellSet<2u> rsqset(am, veclen, contrdepth, 68441551);

      CartesianDerivIterator<2> diter(deriv_order);
      const unsigned int nderiv = diter.range_size();

      const double* A = &(rsqset.R[0][0]);
      const double* B = &(rsqset.R[1][0]);
      LIBINT2_REF_REALTYPE Aref[3];
      for (int i = 0; i < 3; ++i) Aref[i] = A[i];
      LIBINT2_REF_REALTYPE Bref[3];
      for (int i = 0; i < 3; ++i) Bref[i] = B[i];

      prep_libint2(inteval, rsqset, 0, deriv_order);

      cout << "Testing (" << am_to_symbol(am[0]) << "|" << am_to_symbol(am[1])
           << ") ";
      if (deriv_order > 0) {
        cout << " deriv order = " << deriv_order;
      }
      cout << endl;

      double scale_target = 1.0;
#if LIBINT_ACCUM_INTS
      // if accumulating integrals, zero out first, then compute twice
      inteval[0].zero_out_targets = 1;
      scale_target = 0.5;
      if (deriv_order == 0)
        LIBINT2_PREFIXED_NAME(libint2_build_2eri)[am[0]][am[1]](&inteval[0]);
#if INCLUDE_ERI2 >= 1
      if (deriv_order == 1)
        LIBINT2_PREFIXED_NAME(libint2_build_2eri1)[am[0]][am[1]](&inteval[0]);
#endif
#if INCLUDE_ERI2 >= 2
      if (deriv_order == 2)
        LIBINT2_PREFIXED_NAME(libint2_build_2eri2)[am[0]][am[1]](&inteval[0]);
#endif
#if INCLUDE_ERI2 >= 3
      if (deriv_order == 3)
        LIBINT2_PREFIXED_NAME(libint2_build_2eri3)[am[0]][am[1]](&inteval[0]);
#endif
#if INCLUDE_ERI2 >= 4
      if (deriv_order == 4)
        LIBINT2_PREFIXED_NAME(libint2_build_2eri4)[am[0]][am[1]](&inteval[0]);
#endif
#endif
#if LIBINT_CONTRACTED_INTS
      inteval[0].contrdepth = contrdepth2;
#endif
      if (deriv_order == 0)
        LIBINT2_PREFIXED_NAME(libint2_build_2eri)[am[0]][am[1]](&inteval[0]);
#if INCLUDE_ERI2 >= 1
      if (deriv_order == 1)
        LIBINT2_PREFIXED_NAME(libint2_build_2eri1)[am[0]][am[1]](&inteval[0]);
#endif
#if INCLUDE_ERI2 >= 2
      if (deriv_order == 2)
        LIBINT2_PREFIXED_NAME(libint2_build_2eri2)[am[0]][am[1]](&inteval[0]);
#endif
#if INCLUDE_ERI2 >= 3
      if (deriv_order == 3)
        LIBINT2_PREFIXED_NAME(libint2_build_2eri3)[am[0]][am[1]](&inteval[0]);
#endif
#if INCLUDE_ERI2 >= 4
      if (deriv_order == 4)
        LIBINT2_PREFIXED_NAME(libint2_build_2eri4)[am[0]][am[1]](&inteval[0]);
#endif

      bool success = true;
      int ij = 0;

      int l0, m0, n0;
      FOR_CART(l0, m0, n0, am[0])

      int l1, m1, n1;
      FOR_CART(l1, m1, n1, am[1])

      for (uint v = 0; v < veclen; v++) {
        std::vector<LIBINT2_REF_REALTYPE> ref_eri(nderiv, 0.0);

        uint p01 = 0;
        for (uint p0 = 0; p0 < contrdepth; p0++) {
          for (uint p1 = 0; p1 < contrdepth; p1++, p01++) {
            const LIBINT2_REF_REALTYPE alpha0 = rsqset.exp[0][v][p0];
            const LIBINT2_REF_REALTYPE alpha1 = rsqset.exp[1][v][p1];

            const LIBINT2_REF_REALTYPE c0 = rsqset.coef[0][v][p0];
            const LIBINT2_REF_REALTYPE c1 = rsqset.coef[1][v][p1];
            const LIBINT2_REF_REALTYPE c01 = c0 * c1;

            CartesianDerivIterator<2> diter(deriv_order);
            bool last_deriv = false;
            unsigned int di = 0;
            do {
              // convert 2-center deriv indices into 4-center deriv indices
              unsigned int deriv_level[12];
              auto* deriv2_level = &(*diter)[0];
              std::copy(deriv2_level, deriv2_level + 3, deriv_level);
              std::fill(deriv_level + 3, deriv_level + 6, 0u);
              std::copy(deriv2_level + 3, deriv2_level + 6, deriv_level + 6);
              std::fill(deriv_level + 9, deriv_level + 12, 0u);

              ref_eri[di++] +=
                  c01 * eri(deriv_level, l0, m0, n0, alpha0, Aref, 0u, 0u, 0u,
                            0.0, Aref, l1, m1, n1, alpha1, Bref, 0u, 0u, 0u,
                            0.0, Bref, 0);
              last_deriv = diter.last();
              if (!last_deriv) diter.next();
            } while (!last_deriv);
          }
        }

        //
        // extract Libint integrals
        //
        std::vector<LIBINT2_REALTYPE> new_eri;
        for (auto d = 0; d != nderiv; ++d)
          new_eri.push_back(scale_target *
                            inteval[0].targets[d][ij * veclen + v]);

        for (unsigned int di = 0; di < nderiv; ++di) {
          const LIBINT2_REF_REALTYPE abs_error = abs(ref_eri[di] - new_eri[di]);
          const LIBINT2_REF_REALTYPE relabs_error =
              abs(abs_error / ref_eri[di]);
          if (relabs_error > RELATIVE_DEVIATION_THRESHOLD &&
              abs_error > absolute_derivation_threshold(deriv_order)) {
            std::cout << "Elem " << ij << " di= " << di << " v=" << v
                      << " : ref = " << ref_eri[di]
                      << " libint = " << new_eri[di]
                      << " relabs_error = " << relabs_error << endl;
            success = false;
          }
        }

      }  // end of vector loop
      ++ij;
      END_FOR_CART
      END_FOR_CART

      cout << (success ? "ok" : "failed") << endl;

      test_success &= success;
    }
  }

  if (deriv_order == 0)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_2eri)(&inteval[0]);
#if INCLUDE_ERI2 >= 1
  if (deriv_order == 1)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_2eri1)(&inteval[0]);
#endif
#if INCLUDE_ERI2 >= 2
  if (deriv_order == 2)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_2eri2)(&inteval[0]);
#endif
#if INCLUDE_ERI2 >= 3
  if (deriv_order == 3)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_2eri3)(&inteval[0]);
#endif
#if INCLUDE_ERI2 >= 4
  if (deriv_order == 4)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_2eri4)(&inteval[0]);
#endif
  free(inteval);

  return test_success;
}
#endif  // INCLUDE_ERI2
