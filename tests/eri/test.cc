#include <iostream>
#include <cmath>

#include <rr.h>
#include <iter.h>
#include <util.h>
#include <policy_spec.h>
#include <global_macros.h>

#include <libint2.h>
#include <test_eri/eri.h>
#include <test_eri/prep_libint2.h>

using namespace std;
using namespace libint2;

const double RELATIVE_DEVIATION_THRESHOLD = 1.0E-9; // indicate failure if any integral differs in relative sense by more than this

// test 4, 3, and 2-center integrals
#ifdef INCLUDE_ERI
void test_4eri(unsigned int deriv_order);
#endif
#ifdef INCLUDE_ERI3
void test_3eri(unsigned int deriv_order);
#endif
#ifdef INCLUDE_ERI2
void test_2eri(unsigned int deriv_order);
#endif

int main(int argc, char** argv) {
  assert(argc == 1 || argc == 2);
  const unsigned int deriv_order = (argc == 2) ? atoi(argv[1]) : 0u;

  LIBINT2_PREFIXED_NAME(libint2_static_init)();

#ifdef INCLUDE_ERI
  test_4eri(deriv_order);
#endif
#ifdef INCLUDE_ERI3
  test_3eri(deriv_order);
#endif
#ifdef INCLUDE_ERI2
  test_2eri(deriv_order);
#endif

  LIBINT2_PREFIXED_NAME(libint2_static_cleanup)();

  return 0;
}

#ifdef INCLUDE_ERI
void test_4eri(unsigned int deriv_order) {

  if (deriv_order > INCLUDE_ERI) return;

  typedef unsigned int uint;

  const uint veclen = LIBINT2_MAX_VECLEN;
#if LIBINT_CONTRACTED_INTS
  const uint contrdepth = 3;
#else
  const uint contrdepth = 1;
#endif
  const uint contrdepth4 = contrdepth * contrdepth * contrdepth * contrdepth;

  unsigned int lmax;
  if (deriv_order == 0) lmax = LIBINT2_MAX_AM_ERI;
#if INCLUDE_ERI >= 1
  if (deriv_order == 1) lmax = LIBINT2_MAX_AM_ERI1;
#endif
#if INCLUDE_ERI >= 2
  if (deriv_order == 2) lmax = LIBINT2_MAX_AM_ERI2;
#endif
  std::vector<Libint_t> inteval(contrdepth4);
  if (deriv_order == 0) LIBINT2_PREFIXED_NAME(libint2_init_eri)(&inteval[0], lmax, 0);
#if INCLUDE_ERI >= 1
  if (deriv_order == 1) LIBINT2_PREFIXED_NAME(libint2_init_eri1)(&inteval[0], lmax, 0);
#endif
#if INCLUDE_ERI >= 2
  if (deriv_order == 2) LIBINT2_PREFIXED_NAME(libint2_init_eri2)(&inteval[0], lmax, 0);
#endif

  for (unsigned int l0 = 0; l0 <= lmax; ++l0) {
    for (unsigned int l1 = 0; l1 <= lmax; ++l1) {
      for (unsigned int l2 = 0; l2 <= lmax; ++l2) {
        for (unsigned int l3 = 0; l3 <= lmax; ++l3) {

          // can compute this? skip, if not
          if (deriv_order == 0 && LIBINT2_PREFIXED_NAME(libint2_build_eri)[l0][l1][l2][l3] == 0)
            continue;
#if INCLUDE_ERI >= 1
          if (deriv_order == 1 && LIBINT2_PREFIXED_NAME(libint2_build_eri1)[l0][l1][l2][l3] == 0)
            continue;
#endif
#if INCLUDE_ERI >= 2
          if (deriv_order == 2 && LIBINT2_PREFIXED_NAME(libint2_build_eri2)[l0][l1][l2][l3] == 0)
            continue;
#endif

          unsigned int am[4];
          am[0] = l0;
          am[1] = l1;
          am[2] = l2;
          am[3] = l3;
          RandomShellSet<4u> rsqset(am, veclen, contrdepth);

          DerivIndexIterator<4> diter(deriv_order);
          const unsigned int nderiv = diter.range_rank();

          CGShell sh0(am[0]);
          CGShell sh1(am[1]);
          CGShell sh2(am[2]);
          CGShell sh3(am[3]);

          const double* A = &(rsqset.R[0][0]);
          const double* B = &(rsqset.R[1][0]);
          const double* C = &(rsqset.R[2][0]);
          const double* D = &(rsqset.R[3][0]);

          typedef SubIteratorBase<CGShell> iter;
          SafePtr<iter> sh0_iter(new iter(sh0));
          SafePtr<iter> sh1_iter(new iter(sh1));
          SafePtr<iter> sh2_iter(new iter(sh2));
          SafePtr<iter> sh3_iter(new iter(sh3));

          prep_libint2(inteval, rsqset, 0);

          cout << "Testing (" << sh0.label() << sh1.label() << "|"
              << sh2.label() << sh3.label() << ") ";
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
            LIBINT2_PREFIXED_NAME(libint2_build_eri)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);
#if INCLUDE_ERI >= 1
          else if (deriv_order == 1)
            LIBINT2_PREFIXED_NAME(libint2_build_eri1)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);
#endif
#if INCLUDE_ERI >= 2
          else if (deriv_order == 2)
            LIBINT2_PREFIXED_NAME(libint2_build_eri2)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);
#endif
#endif
#if LIBINT_CONTRACTED_INTS
          inteval[0].contrdepth = contrdepth4;
#endif
          if (deriv_order == 0)
            LIBINT2_PREFIXED_NAME(libint2_build_eri)[am[0]][am[1]][am[2]][am[3]](
                &inteval[0]);
#if INCLUDE_ERI >= 1
          else if (deriv_order == 1)
            LIBINT2_PREFIXED_NAME(libint2_build_eri1)[am[0]][am[1]][am[2]][am[3]](
                &inteval[0]);
#endif
#if INCLUDE_ERI >= 2
          else if (deriv_order == 2)
            LIBINT2_PREFIXED_NAME(libint2_build_eri2)[am[0]][am[1]][am[2]][am[3]](
                &inteval[0]);
#endif

          bool success = true;
          int ijkl = 0;
          for (sh0_iter->init(); int(*sh0_iter); ++(*sh0_iter)) {
            for (sh1_iter->init(); int(*sh1_iter); ++(*sh1_iter)) {
              for (sh2_iter->init(); int(*sh2_iter); ++(*sh2_iter)) {
                for (sh3_iter->init(); int(*sh3_iter); ++(*sh3_iter), ijkl++) {

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

                  for (uint v = 0; v < veclen; v++) {

                    std::vector<double> ref_eri(nderiv, 0.0);

                    uint p0123 = 0;
                    for (uint p0 = 0; p0 < contrdepth; p0++) {
                      for (uint p1 = 0; p1 < contrdepth; p1++) {
                        for (uint p2 = 0; p2 < contrdepth; p2++) {
                          for (uint p3 = 0; p3 < contrdepth; p3++, p0123++) {

                            const double alpha0 = rsqset.exp[0][v][p0];
                            const double alpha1 = rsqset.exp[1][v][p1];
                            const double alpha2 = rsqset.exp[2][v][p2];
                            const double alpha3 = rsqset.exp[3][v][p3];

                            const double c0 = rsqset.coef[0][v][p0];
                            const double c1 = rsqset.coef[1][v][p1];
                            const double c2 = rsqset.coef[2][v][p2];
                            const double c3 = rsqset.coef[3][v][p3];
                            const double c0123 = c0 * c1 * c2 * c3;

                            DerivIndexIterator<4> diter(deriv_order);
                            bool last_deriv = false;
                            unsigned int di = 0;
                            do {
                              ref_eri[di++] += c0123
                                  * eri(diter.values(), l0, m0, n0, alpha0, A,
                                        l1, m1, n1, alpha1, B, l2, m2, n2,
                                        alpha2, C, l3, m3, n3, alpha3, D, 0);
                              last_deriv = diter.last();
                              if (!last_deriv)
                                diter.next();
                            } while (!last_deriv);

                          }
                        }
                      }
                    }

                    std::vector<double> new_eri;
                    if (deriv_order == 0)
                      new_eri.push_back( scale_target * inteval[0].targets[0][ijkl * veclen + v] );

                    // derivatives w.r.t. center 2 skipped and must be reconstructed using the translational invariance
                    if (deriv_order == 1) {
                      for (unsigned int di = 0; di < nderiv; ++di) {
                        double value;
                        if (di<6)
                          value = scale_target * inteval[0].targets[di][ijkl * veclen + v];
                        else if (di>=9)
                          value = scale_target * inteval[0].targets[di-3][ijkl * veclen + v];
                        else // (di>=6 || di<=8)
                          value = -scale_target * (inteval[0].targets[di-6][ijkl * veclen + v] +
                              inteval[0].targets[di-3][ijkl * veclen + v] +
                              inteval[0].targets[di][ijkl * veclen + v]);
                        new_eri.push_back(value);
                      }
                    }

                    // derivatives w.r.t. center 2 are skipped and must be reconstructed using the translational invariance
                    if (deriv_order == 2) {
                      const unsigned int NumCenters = 4;
                      double libint_eri[3*(NumCenters-1)][3*(NumCenters-1)];
                      for(int di=0,dij=0; di<3*(NumCenters-1); di++) {
                        for(int dj=di; dj<3*(NumCenters-1); dj++, ++dij){
                          const double value = scale_target * inteval[0].targets[dij][ijkl * veclen + v];
                          libint_eri[di][dj] = value;
                          libint_eri[dj][di] = value;
                        }
                      }

                      for(int di=0; di<3*NumCenters; di++) {
                        for(int dj=di; dj<3*NumCenters; dj++){

                          const int ci = di/3;
                          const int cj = dj/3;
                          const int xyzi = di%3;
                          const int xyzj = dj%3;

                          double value = 0.0;

                          // d2/dCidCj = d2/dAidAj + d2/dAidBj + d2/dAidDj + d2/dBidAj + d2/dBidBj + d2/dBidDj + d2/dDidAj + d2/dDidBj + d2/dDidDj
                          if (ci == 2 && cj == 2) {
                            value = libint_eri[xyzi][xyzj] + libint_eri[xyzi][xyzj+3] + libint_eri[xyzi][xyzj+6] +
                                    libint_eri[xyzi+3][xyzj] + libint_eri[xyzi+3][xyzj+3] + libint_eri[xyzi+3][xyzj+6] +
                                    libint_eri[xyzi+6][xyzj] + libint_eri[xyzi+6][xyzj+3] + libint_eri[xyzi+6][xyzj+6];
                          }
                          // d2/dCidDj = - d2/dAidDj - d2/dBidDj - -d2/dDidDj
                          else if (ci == 2) {
                            value = -libint_eri[xyzi][dj-3] - libint_eri[xyzi+3][dj-3] - libint_eri[xyzi+6][dj-3] ;
                          }
                          // d2/dXidCj = - d2/dXidAj - d2/dXidBj - -d2/dXidDj (X=A,B)
                          else if (cj == 2) {
                            value = -libint_eri[di][xyzj] - libint_eri[di][xyzj+3] - libint_eri[di][xyzj+6] ;
                          }
                          // d2/dDidDj
                          else if (ci == 3) {
                            value = libint_eri[di-3][dj-3];
                          }
                          // d2/dXidDj (X=A,B)
                          else if (cj == 3) {
                            value = libint_eri[di][dj-3];
                          }
                          // d2/dXidYj (X=A,B)
                          else {
                            value = libint_eri[di][dj];
                          }

                          new_eri.push_back(value);
                        }
                      }
                    }

                    for (unsigned int di = 0; di < nderiv; ++di) {
                      if (fabs((ref_eri[di] - new_eri[di]) / new_eri[di]) > RELATIVE_DEVIATION_THRESHOLD) {
                        std::cout << "Elem " << ijkl << " di= " << di << " v="
                            << v << " : eri.cc = " << ref_eri[di]
                            << " libint = " << new_eri[di] << endl;
                        success = false;
                      }
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

  if (deriv_order == 0)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);
#if INCLUDE_ERI >= 1
  if (deriv_order == 1)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);
#endif
#if INCLUDE_ERI >= 2
  if (deriv_order == 2)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);
#endif
}
#endif // INCLUDE_ERI

#ifdef INCLUDE_ERI3
void test_3eri(unsigned int deriv_order) {

  if (deriv_order > INCLUDE_ERI3) return;

  typedef unsigned int uint;

  const uint veclen = LIBINT2_MAX_VECLEN;
#if LIBINT_CONTRACTED_INTS
  const uint contrdepth = 3;
#else
  const uint contrdepth = 1;
#endif
  const uint contrdepth3 = contrdepth * contrdepth * contrdepth;

  unsigned int lmax;
  if (deriv_order == 0) lmax = LIBINT2_MAX_AM_3ERI;
#if INCLUDE_ERI3 >= 1
  if (deriv_order == 1) lmax = LIBINT2_MAX_AM_3ERI1;
#endif
#if INCLUDE_ERI3 >= 2
  if (deriv_order == 2) lmax = LIBINT2_MAX_AM_3ERI2;
#endif
  std::vector<Libint_t> inteval(contrdepth3);
  if (deriv_order == 0)
    LIBINT2_PREFIXED_NAME(libint2_init_3eri)(&inteval[0], lmax, 0);
#if INCLUDE_ERI3 >= 1
  if (deriv_order == 1)
    LIBINT2_PREFIXED_NAME(libint2_init_3eri1)(&inteval[0], lmax, 0);
#endif
#if INCLUDE_ERI3 >= 2
  if (deriv_order == 2)
    LIBINT2_PREFIXED_NAME(libint2_init_3eri2)(&inteval[0], lmax, 0);
#endif

  for (unsigned int l0 = 0; l0 <= lmax; ++l0) {
    for (unsigned int l1 = 0; l1 <= lmax; ++l1) {
      for (unsigned int l2 = 0; l2 <= lmax; ++l2) {

        // can compute this? skip, if not
        if (deriv_order == 0 && LIBINT2_PREFIXED_NAME(libint2_build_3eri)[l0][l1][l2] == 0)
          continue;
#if INCLUDE_ERI3 >= 1
        if (deriv_order == 1 && LIBINT2_PREFIXED_NAME(libint2_build_3eri1)[l0][l1][l2] == 0)
          continue;
#endif
#if INCLUDE_ERI3 >= 2
        if (deriv_order == 2 && LIBINT2_PREFIXED_NAME(libint2_build_3eri2)[l0][l1][l2] == 0)
          continue;
#endif

        unsigned int am[3];
        am[0] = l0;
        am[1] = l1;
        am[2] = l2;
        RandomShellSet<3u> rsqset(am, veclen, contrdepth);

        DerivIndexIterator<3> diter(deriv_order);
        const unsigned int nderiv = diter.range_rank();

        CGShell sh0(am[0]);
        CGShell sh1(am[1]);
        CGShell sh2(am[2]);

        const double* A = &(rsqset.R[0][0]);
        const double* B = &(rsqset.R[1][0]);
        const double* C = &(rsqset.R[2][0]);

        typedef SubIteratorBase<CGShell> iter;
        SafePtr<iter> sh0_iter(new iter(sh0));
        SafePtr<iter> sh1_iter(new iter(sh1));
        SafePtr<iter> sh2_iter(new iter(sh2));

        prep_libint2(inteval, rsqset, 0);

        cout << "Testing (" << sh0.label() << "|" << sh1.label() << sh2.label()
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
          LIBINT2_PREFIXED_NAME(libint2_build_3eri)[am[0]][am[1]][am[2]](&inteval[0]);
#if INCLUDE_ERI3 >= 1
        if (deriv_order == 1)
          LIBINT2_PREFIXED_NAME(libint2_build_3eri1)[am[0]][am[1]][am[2]](&inteval[0]);
#endif
#if INCLUDE_ERI3 >= 2
        if (deriv_order == 2)
          LIBINT2_PREFIXED_NAME(libint2_build_3eri2)[am[0]][am[1]][am[2]](&inteval[0]);
#endif
#endif
#if LIBINT_CONTRACTED_INTS
        inteval[0].contrdepth = contrdepth3;
#endif
        if (deriv_order == 0)
          LIBINT2_PREFIXED_NAME(libint2_build_3eri)[am[0]][am[1]][am[2]](&inteval[0]);
#if INCLUDE_ERI3 >= 1
        if (deriv_order == 1)
          LIBINT2_PREFIXED_NAME(libint2_build_3eri1)[am[0]][am[1]][am[2]](&inteval[0]);
#endif
#if INCLUDE_ERI3 >= 2
        if (deriv_order == 2)
          LIBINT2_PREFIXED_NAME(libint2_build_3eri2)[am[0]][am[1]][am[2]](&inteval[0]);
#endif

        bool success = true;
        int ijk = 0;
        for (sh0_iter->init(); int(*sh0_iter); ++(*sh0_iter)) {
          for (sh1_iter->init(); int(*sh1_iter); ++(*sh1_iter)) {
            for (sh2_iter->init(); int(*sh2_iter); ++(*sh2_iter), ijk++) {

              CGF bf0 = sh0_iter->elem();
              CGF bf1 = sh1_iter->elem();
              CGF bf2 = sh2_iter->elem();

              uint l0 = bf0.qn(0);
              uint m0 = bf0.qn(1);
              uint n0 = bf0.qn(2);
              uint l1 = bf1.qn(0);
              uint m1 = bf1.qn(1);
              uint n1 = bf1.qn(2);
              uint l2 = bf2.qn(0);
              uint m2 = bf2.qn(1);
              uint n2 = bf2.qn(2);

              for (uint v = 0; v < veclen; v++) {

                std::vector<double> ref_eri(nderiv, 0.0);

                uint p012 = 0;
                for (uint p0 = 0; p0 < contrdepth; p0++) {
                  for (uint p1 = 0; p1 < contrdepth; p1++) {
                    for (uint p2 = 0; p2 < contrdepth; p2++, p012++) {

                      const double alpha0 = rsqset.exp[0][v][p0];
                      const double alpha1 = rsqset.exp[1][v][p1];
                      const double alpha2 = rsqset.exp[2][v][p2];

                      const double c0 = rsqset.coef[0][v][p0];
                      const double c1 = rsqset.coef[1][v][p1];
                      const double c2 = rsqset.coef[2][v][p2];
                      const double c012 = c0 * c1 * c2;

                      DerivIndexIterator<3> diter(deriv_order);
                      bool last_deriv = false;
                      unsigned int di = 0;
                      do {
                        // convert 3-center deriv indices into 4-center deriv indices
                        unsigned int deriv_level[12];
                        std::copy(diter.values(), diter.values()+3, deriv_level);
                        std::fill(deriv_level+3, deriv_level+6, 0u);
                        std::copy(diter.values()+3, diter.values()+9, deriv_level+6);

                        ref_eri[di++] += c012
                            * eri(deriv_level, l0, m0, n0, alpha0, A, 0u, 0u,
                                  0u, 0.0, A, l1, m1, n1, alpha1, B, l2, m2, n2,
                                  alpha2, C, 0);
                        last_deriv = diter.last();
                        if (!last_deriv)
                          diter.next();
                      } while (!last_deriv);

                    }
                  }
                }

                std::vector<double> new_eri;

                if (deriv_order == 0)
                  new_eri.push_back( scale_target * inteval[0].targets[0][ijk * veclen + v] );

                // derivatives w.r.t. center 0 are skipped and must be reconstructed using the translational invariance
                if (deriv_order == 1) {
                  for (unsigned int di = 0; di < nderiv; ++di) {
                    if (di>=3)
                      new_eri.push_back( scale_target * inteval[0].targets[di-3][ijk * veclen + v] );
                    else // (di<3)
                      new_eri.push_back( -scale_target * (inteval[0].targets[di][ijk * veclen + v] +
                                                          inteval[0].targets[di+3][ijk * veclen + v])
                                                          );
                  }
                }

                // derivatives w.r.t. center 0 are skipped and must be reconstructed using the translational invariance
                if (deriv_order == 2) {
                  const unsigned int NumCenters = 3;
                  double libint_eri[3*(NumCenters-1)][3*(NumCenters-1)];
                  for(int di=0,dij=0; di<3*(NumCenters-1); di++) {
                    for(int dj=di; dj<3*(NumCenters-1); dj++, ++dij){
                      const double value = scale_target * inteval[0].targets[dij][ijk * veclen + v];
                      libint_eri[di][dj] = value;
                      libint_eri[dj][di] = value;
                    }
                  }

                  for(int di=0; di<3*NumCenters; di++) {
                    for(int dj=di; dj<3*NumCenters; dj++){

                      const int ci = di/3;
                      const int cj = dj/3;
                      const int xyzi = di%3;
                      const int xyzj = dj%3;

                      double value = 0.0;

                      // d2/dAidAj = d2/dBidBj + d2/dBidCj + d2/dCidBj + d2/dCidCj
                      if (ci == 0 && cj == 0) {
                        value = libint_eri[xyzi][xyzj] + libint_eri[xyzi][3+xyzj] + libint_eri[xyzi+3][xyzj] + libint_eri[3+xyzi][3+xyzj];
                      }
                      // d2/dAidXj = - d2/dBidXj - -d2/dCidXj
                      else if (ci == 0) {
                        value = -libint_eri[xyzi][dj-3] - libint_eri[3+xyzi][dj-3];
                      }
                      // rest
                      else {
                        value = libint_eri[di-3][dj-3];
                      }

                      new_eri.push_back(value);
                    }
                  }
                }

                for (unsigned int di = 0; di < nderiv; ++di) {
                  if (fabs((ref_eri[di] - new_eri[di]) / new_eri[di]) > RELATIVE_DEVIATION_THRESHOLD) {
                    std::cout << "Elem " << ijk << " di= " << di << " v=" << v
                        << " : eri.cc = " << ref_eri[di] << " libint = "
                        << new_eri[di] << endl;
                    success = false;
                  }
                }

              } // end of vector loop
            }
          }
        }

        cout << (success ? "ok" : "failed") << endl;

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
}
#endif // INCLUDE_ERI3

#ifdef INCLUDE_ERI2
void test_2eri(unsigned int deriv_order) {

  if (deriv_order > INCLUDE_ERI2) return;

  typedef unsigned int uint;

  const uint veclen = LIBINT2_MAX_VECLEN;
#if LIBINT_CONTRACTED_INTS
  const uint contrdepth = 3;
#else
  const uint contrdepth = 1;
#endif
  const uint contrdepth2 = contrdepth * contrdepth;

  unsigned int lmax;
  if (deriv_order == 0) lmax = LIBINT2_MAX_AM_2ERI;
#if INCLUDE_ERI2 >= 1
  if (deriv_order == 1) lmax = LIBINT2_MAX_AM_2ERI1;
#endif
#if INCLUDE_ERI2 >= 2
  if (deriv_order == 2) lmax = LIBINT2_MAX_AM_2ERI2;
#endif
  std::vector<Libint_t> inteval(contrdepth2);
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

  for (unsigned int l0 = 0; l0 <= lmax; ++l0) {
    for (unsigned int l1 = 0; l1 <= lmax; ++l1) {

        // can compute this? skip, if not
        if (deriv_order == 0 && LIBINT2_PREFIXED_NAME(libint2_build_2eri)[l0][l1] == 0)
          continue;
#if INCLUDE_ERI2 >= 1
        if (deriv_order == 1 && LIBINT2_PREFIXED_NAME(libint2_build_2eri1)[l0][l1] == 0)
          continue;
#endif
#if INCLUDE_ERI2 >= 2
        if (deriv_order == 2 && LIBINT2_PREFIXED_NAME(libint2_build_2eri2)[l0][l1] == 0)
          continue;
#endif

        unsigned int am[2];
        am[0] = l0;
        am[1] = l1;
        RandomShellSet<2u> rsqset(am, veclen, contrdepth);

        DerivIndexIterator<2> diter(deriv_order);
        const unsigned int nderiv = diter.range_rank();

        CGShell sh0(am[0]);
        CGShell sh1(am[1]);

        const double* A = &(rsqset.R[0][0]);
        const double* B = &(rsqset.R[1][0]);

        typedef SubIteratorBase<CGShell> iter;
        SafePtr<iter> sh0_iter(new iter(sh0));
        SafePtr<iter> sh1_iter(new iter(sh1));

        prep_libint2(inteval, rsqset, 0);

        cout << "Testing (" << sh0.label() << "|" << sh1.label() << ") ";
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

        bool success = true;
        int ij = 0;
        for (sh0_iter->init(); int(*sh0_iter); ++(*sh0_iter)) {
          for (sh1_iter->init(); int(*sh1_iter); ++(*sh1_iter), ij++) {

              CGF bf0 = sh0_iter->elem();
              CGF bf1 = sh1_iter->elem();

              uint l0 = bf0.qn(0);
              uint m0 = bf0.qn(1);
              uint n0 = bf0.qn(2);
              uint l1 = bf1.qn(0);
              uint m1 = bf1.qn(1);
              uint n1 = bf1.qn(2);

              for (uint v = 0; v < veclen; v++) {

                std::vector<double> ref_eri(nderiv, 0.0);

                uint p01 = 0;
                for (uint p0 = 0; p0 < contrdepth; p0++) {
                  for (uint p1 = 0; p1 < contrdepth; p1++, p01++) {

                      const double alpha0 = rsqset.exp[0][v][p0];
                      const double alpha1 = rsqset.exp[1][v][p1];

                      const double c0 = rsqset.coef[0][v][p0];
                      const double c1 = rsqset.coef[1][v][p1];
                      const double c01 = c0 * c1;

                      DerivIndexIterator<2> diter(deriv_order);
                      bool last_deriv = false;
                      unsigned int di = 0;
                      do {
                        // convert 2-center deriv indices into 4-center deriv indices
                        unsigned int deriv_level[12];
                        std::copy(diter.values(), diter.values()+3, deriv_level);
                        std::fill(deriv_level+3, deriv_level+6, 0u);
                        std::copy(diter.values()+3, diter.values()+6, deriv_level+6);
                        std::fill(deriv_level+9, deriv_level+12, 0u);

                        ref_eri[di++] += c01
                            * eri(deriv_level, l0, m0, n0, alpha0, A, 0u, 0u,
                                  0u, 0.0, A, l1, m1, n1, alpha1, B, 0u, 0u,
                                  0u, 0.0, B, 0);
                        last_deriv = diter.last();
                        if (!last_deriv)
                          diter.next();
                      } while (!last_deriv);

                  }
                }

                std::vector<double> new_eri;

                if (deriv_order == 0)
                  new_eri.push_back(scale_target * inteval[0].targets[0][ij * veclen + v]);

                // derivatives w.r.t. center 0 are skipped and must be reconstructed using the translational invariance
                if (deriv_order == 1) {
                  for (unsigned int di = 0; di < nderiv; ++di) {
                    if (di>=3)
                      new_eri.push_back(scale_target * inteval[0].targets[di-3][ij * veclen + v]);
                    else // (di<3)
                      new_eri.push_back(-scale_target * inteval[0].targets[di][ij * veclen + v]);
                  }
                }

                // derivatives w.r.t. center 0 are skipped and must be reconstructed using the translational invariance
                if (deriv_order == 2) {

                  double libint_eri[3][3];
                  for(int di=0,dij=0; di<3; di++) {
                    for(int dj=di; dj<3; dj++, ++dij){
                      const double value = scale_target * inteval[0].targets[dij][ij * veclen + v];
                      libint_eri[di][dj] = value;
                      libint_eri[dj][di] = value;
                    }
                  }

                  for(int di=0; di<6; di++) {
                    for(int dj=di; dj<6; dj++){

                      const int ci = di/3;
                      const int cj = dj/3;
                      const int xyzi = di%3;
                      const int xyzj = dj%3;

                      double value = 0.0;

                      // d2/dAidAj
                      if (ci == 0 && cj == 0) {
                        value = libint_eri[xyzi][xyzj];
                      }
                      // d2/dAidBj
                      else if (ci == 0) {
                        value = -libint_eri[xyzi][xyzj];
                      }
                      // rest
                      else {
                        value = libint_eri[xyzi][xyzj];
                      }

                      new_eri.push_back(value);
                    }
                  }
                }

                for (unsigned int di = 0; di < nderiv; ++di) {
                  if (fabs((ref_eri[di] - new_eri[di]) / new_eri[di]) > RELATIVE_DEVIATION_THRESHOLD) {
                    std::cout << "Elem " << ij << " di= " << di << " v=" << v
                        << " : eri.cc = " << ref_eri[di] << " libint = "
                        << new_eri[di] << endl;
                    success = false;
                  }
                }

              } // end of vector loop
          }
        }

        cout << (success ? "ok" : "failed") << endl;

    }
  }

  if (deriv_order == 0) LIBINT2_PREFIXED_NAME(libint2_cleanup_2eri)(&inteval[0]);
#if INCLUDE_ERI2 >= 1
  if (deriv_order == 1) LIBINT2_PREFIXED_NAME(libint2_cleanup_2eri1)(&inteval[0]);
#endif
#if INCLUDE_ERI2 >= 2
  if (deriv_order == 2) LIBINT2_PREFIXED_NAME(libint2_cleanup_2eri2)(&inteval[0]);
#endif

}
#endif // INCLUDE_ERI2
