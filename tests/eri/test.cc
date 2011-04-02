#include <iostream>
#include <cmath>

#include <rr.h>
#include <iter.h>
#include <policy_spec.h>
#include <global_macros.h>

#include <libint2.h>
#include <test_eri/eri.h>
#include <test_eri/prep_libint2.h>

using namespace std;
using namespace libint2;

int main(int argc, char** argv) {

  typedef unsigned int uint;

  LIBINT2_PREFIXED_NAME(libint2_static_init)();

  const uint veclen = LIBINT2_MAX_VECLEN;
  double alpha[4] = { 0.5, 1.0, 1.5, 2.0 };
  double A[3] = { 1.0, 2.0, 3.0 };
  double B[3] = { 1.5, 2.5, 3.5 };
  double C[3] = { 4.0, 2.0, 0.0 };
  double D[3] = { 3.0, 3.0, 1.0 };

  const double ratio = 1.5;
  std::vector<double> alpha1;
  std::vector<double> alpha2;
  std::vector<double> alpha3;
  std::vector<double> alpha4;
  for (uint v = 0; v < veclen; v++) {
    const double scale = pow(ratio, static_cast<double> (v));
    alpha1.push_back(alpha[0] * scale);
    alpha2.push_back(alpha[1] * scale);
    alpha3.push_back(alpha[2] * scale);
    alpha4.push_back(alpha[3] * scale);
  }

  const unsigned int lmax = LIBINT2_MAX_AM_ERI;
  Libint_t inteval;
  LIBINT2_PREFIXED_NAME(libint2_init_eri)(&inteval, lmax, 0);

  for (unsigned int l0 = 0; l0 <= lmax; ++l0) {
    for (unsigned int l1 = 0; l1 <= lmax; ++l1) {
      for (unsigned int l2 = 0; l2 <= lmax; ++l2) {
        for (unsigned int l3 = 0; l3 <= lmax; ++l3) {

          // can compute this? skip, if not
          if (LIBINT2_PREFIXED_NAME(libint2_build_eri)[l0][l1][l2][l3] == 0)
            continue;

          unsigned int am[4];
          am[0] = l0;
          am[1] = l1;
          am[2] = l2;
          am[3] = l3;

#if USE_BRAKET_H
          CGShell sh0(&(am[0]));
          CGShell sh1(&(am[1]));
          CGShell sh2(&(am[2]));
          CGShell sh3(&(am[3]));
#else
          SafePtr<CGShell> sh0(new CGShell(&(am[0])));
          SafePtr<CGShell> sh1(new CGShell(&(am[1])));
          SafePtr<CGShell> sh2(new CGShell(&(am[2])));
          SafePtr<CGShell> sh3(new CGShell(&(am[3])));
#endif

          typedef SubIteratorBase<CGShell> iter;
          SafePtr<iter> sh0_iter(new iter(sh0));
          SafePtr<iter> sh1_iter(new iter(sh1));
          SafePtr<iter> sh2_iter(new iter(sh2));
          SafePtr<iter> sh3_iter(new iter(sh3));

          prep_libint2(&inteval, am[0], alpha1, A, am[1], alpha2, B, am[2],
                       alpha3, C, am[3], alpha4, D, 0, veclen);

#if USE_BRAKET_H
          cout << "Testing (" << sh0.label() << sh1.label() << "|"
              << sh2.label() << sh3.label() << ") ... ";
#else
          cout << "Testing (" << sh0->label() << sh1->label()
          << "|" << sh2->label() << sh3->label() << ") ... ";
#endif

          double scale_target = 1.0;
#if LIBINT_ACCUM_INTS
          // if accumulating integrals, zero out first, then compute twice
          erieval->zero_out_targets = 1;
          scale_target = 0.5;
          LIBINT2_PREFIXED_NAME(libint2_build_eri)[am[0]][am[1]][am[2]][am[3]](&inteval);
#endif
          LIBINT2_PREFIXED_NAME(libint2_build_eri)[am[0]][am[1]][am[2]][am[3]](&inteval);

          bool success = true;
          int ijkl = 0;
          for (sh0_iter->init(); int(*sh0_iter); ++(*sh0_iter)) {
            for (sh1_iter->init(); int(*sh1_iter); ++(*sh1_iter)) {
              for (sh2_iter->init(); int(*sh2_iter); ++(*sh2_iter)) {
                for (sh3_iter->init(); int(*sh3_iter); ++(*sh3_iter), ijkl++) {

#if USE_BRAKET_H
                  CGF bf0 = sh0_iter->elem();
                  CGF bf1 = sh1_iter->elem();
                  CGF bf2 = sh2_iter->elem();
                  CGF bf3 = sh3_iter->elem();

                  uint l0 = bf0.qn(0);
                  uint m0 = bf0.qn(1);
                  uint n0 = bf0.qn(2);
                  uint l1 = bf1.qn(0);
                  uint m1 = bf1.qn(1);
                  uint n1 = bf1.qn(2);
                  uint l2 = bf2.qn(0);
                  uint m2 = bf2.qn(1);
                  uint n2 = bf2.qn(2);
                  uint l3 = bf3.qn(0);
                  uint m3 = bf3.qn(1);
                  uint n3 = bf3.qn(2);
#else
                  SafePtr<CGF> bf0 = sh0_iter->elem();
                  SafePtr<CGF> bf1 = sh1_iter->elem();
                  SafePtr<CGF> bf2 = sh2_iter->elem();
                  SafePtr<CGF> bf3 = sh3_iter->elem();

                  uint l0 = bf0->qn(0);
                  uint m0 = bf0->qn(1);
                  uint n0 = bf0->qn(2);
                  uint l1 = bf1->qn(0);
                  uint m1 = bf1->qn(1);
                  uint n1 = bf1->qn(2);
                  uint l2 = bf2->qn(0);
                  uint m2 = bf2->qn(1);
                  uint n2 = bf2->qn(2);
                  uint l3 = bf3->qn(0);
                  uint m3 = bf3->qn(1);
                  uint n3 = bf3->qn(2);
#endif

                  for (uint v = 0; v < veclen; v++) {

                    double ref_eri = eri(l0, m0, n0, alpha1[v], A, l1, m1, n1,
                                         alpha2[v], B, l2, m2, n2, alpha3[v],
                                         C, l3, m3, n3, alpha4[v], D, 0);

                    double new_eri = scale_target * inteval.targets[0][ijkl
                        * veclen + v];

                    if (fabs((ref_eri - new_eri) / new_eri) > 1.0E-8) {
#if 0
                      std::cout << std::endl << "Elem " << ijkl << " v=" << v
                          << " : eri.cc = " << ref_eri << " libint = "
                          << new_eri;
#endif
                      success = false;
                    }
                  } // end of vector loop
                }
              }
            }
          }

          cout << (success ? "ok" : "failed") << endl;

        }
      }
    }
  }

  LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval);

  return 0;
}

