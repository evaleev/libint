#include <libint2.h>
#include <math.h>
#include <stdio.h>

unsigned int maxfac = 50;
double df[150];
double F[12];

void calc_f(double* F, int n, double t)
{
  int i, m, k;
  int m2;
  double t2;
  double num;
  double sum;
  double term1, term2;
  double et;
  double K = 1.0/M_2_SQRTPI;
  unsigned int maxfac = 100;

  if (t>20.0){   /* For big t's do upward recursion */
    t2 = 2*t;
    et = exp(-t);
    t = sqrt(t);
    F[0] = K*erf(t)/t;
    for(m=0; m<=n-1; m++){
      F[m+1] = ((2*m + 1)*F[m] - et)/(t2);
    }
  }
  else {        /* For smaller t's compute F with highest n using
                   asymptotic series (see I. Shavitt in
                   Methods in Computational Physics, ed. B. Alder eta l,
                   vol 2, 1963, page 8) */
    et = exp(-t);
    t2 = 2*t;
    m2 = 2*n;
    num = df[m2];
    i=0;
    sum = 1.0/(m2+1);
    do{
      i++;
      num = num*t2;
      term1 = num/df[m2+2*i+2];
      sum += term1;
    } while (term1 > 1e-12 && i < maxfac);
    F[n] = sum*et;
    for(m=n-1;m>=0;m--){        /* And then do downward recursion */
      F[m] = (t2*F[m+1] + et)/(2*m+1);
    }
  }
}


/** This function evaluates ERI over 4 primitive Gaussian shells.
    See tests/eri/test.cc for an example of how to deal with
    contracted Gaussians.

    For simplicity, many details are omitted here, e.g. normalization.
  */
void
compute_eri(Libint_t* erieval,
            unsigned int am1, double alpha1, double* A,
            unsigned int am2, double alpha2, double* B,
            unsigned int am3, double alpha3, double* C,
            unsigned int am4, double alpha4, double* D
           )
{
  /* I will assume that libint2_static_init() and libint2_init_eri(&erieval,max_am,0) had been called elsewhere! */

  //-------------------------------------------------------------------------------------

  double gammap, Px, Py, Pz, PAx, PAy, PAz, PBx, PBy, PBz, AB2;
  double gammaq, Qx, Qy, Qz, QCx, QCy, QCz, QDx, QDy, QDz, CD2;
  double gammapq, PQx, PQy, PQz, PQ2, Wx, Wy, Wz;
  double K1, K2, pfac;
  unsigned int am;
  double* eri_shell_set;
  unsigned int n1, n2, n3, n4;
  int a, b, c, d;

  //
  // Compute requisite data -- many of these quantities would be precomputed
  // for all nonnegligible shell pairs somewhere else
  //
  gammap = alpha1 + alpha2;
  Px = (alpha1*A[0] + alpha2*B[0])/gammap;
  Py = (alpha1*A[1] + alpha2*B[1])/gammap;
  Pz = (alpha1*A[2] + alpha2*B[2])/gammap;
  PAx = Px - A[0];
  PAy = Py - A[1];
  PAz = Pz - A[2];
  PBx = Px - B[0];
  PBy = Py - B[1];
  PBz = Pz - B[2];
  AB2 = (A[0]-B[0])*(A[0]-B[0])
      + (A[1]-B[1])*(A[1]-B[1])
      + (A[2]-B[2])*(A[2]-B[2]);

  erieval->PA_x[0] = PAx;
  erieval->PA_y[0] = PAy;
  erieval->PA_z[0] = PAz;
  erieval->AB_x[0] = A[0] - B[0];
  erieval->AB_y[0] = A[1] - B[1];
  erieval->AB_z[0] = A[2] - B[2];
  erieval->oo2z[0] = 0.5/gammap;

  gammaq = alpha3 + alpha4;
  gammapq = gammap*gammaq/(gammap+gammaq);
  Qx = (alpha3*C[0] + alpha4*D[0])/gammaq;
  Qy = (alpha3*C[1] + alpha4*D[1])/gammaq;
  Qz = (alpha3*C[2] + alpha4*D[2])/gammaq;
  QCx = Qx - C[0];
  QCy = Qy - C[1];
  QCz = Qz - C[2];
  QDx = Qx - D[0];
  QDy = Qy - D[1];
  QDz = Qz - D[2];
  CD2 = (C[0]-D[0])*(C[0]-D[0])
      + (C[1]-D[1])*(C[1]-D[1])
      + (C[2]-D[2])*(C[2]-D[2]);

  erieval->QC_x[0] = QCx;
  erieval->QC_y[0] = QCy;
  erieval->QC_z[0] = QCz;
  erieval->CD_x[0] = C[0] - D[0];
  erieval->CD_y[0] = C[1] - D[1];
  erieval->CD_z[0] = C[2] - D[2];
  erieval->oo2e[0] = 0.5/gammaq;

  PQx = Px - Qx;
  PQy = Py - Qy;
  PQz = Pz - Qz;
  PQ2 = PQx*PQx + PQy*PQy + PQz*PQz;
  Wx = (gammap*Px + gammaq*Qx)/(gammap+gammaq);
  Wy = (gammap*Py + gammaq*Qy)/(gammap+gammaq);
  Wz = (gammap*Pz + gammaq*Qz)/(gammap+gammaq);

  erieval->WP_x[0] = Wx - Px;
  erieval->WP_y[0] = Wy - Py;
  erieval->WP_z[0] = Wz - Pz;
  erieval->WQ_x[0] = Wx - Qx;
  erieval->WQ_y[0] = Wy - Qy;
  erieval->WQ_z[0] = Wz - Qz;
  erieval->oo2ze[0] = 0.5/(gammap+gammaq);
  erieval->roz[0] = gammapq/gammap;
  erieval->roe[0] = gammapq/gammaq;

  K1 = exp(-alpha1*alpha2*AB2/gammap);
  K2 = exp(-alpha3*alpha4*CD2/gammaq);
  pfac = 2*pow(M_PI,2.5)*K1*K2/(gammap*gammaq*sqrt(gammap+gammaq));

  //
  // evaluate Boys function F_m for all m in [0,am]
  //
  am = am1 + am2 + am3 + am4;
  calc_f(F, PQ2*gammapq, am);

  // (00|00)^m = pfac * F_m
  erieval->LIBINT_T_SS_EREP_SS(0)[0] = pfac*F[0];
  erieval->LIBINT_T_SS_EREP_SS(1)[0] = pfac*F[1];
  erieval->LIBINT_T_SS_EREP_SS(2)[0] = pfac*F[2];
  erieval->LIBINT_T_SS_EREP_SS(3)[0] = pfac*F[3];
  erieval->LIBINT_T_SS_EREP_SS(4)[0] = pfac*F[4];

  // compute ERIs
  libint2_build_eri[am1][am2][am3][am4](erieval);

  // Print out the integrals
  eri_shell_set = erieval->targets[0];
  n1 = (am1 + 1) * (am1 + 2)/2;
  n2 = (am2 + 1) * (am2 + 2)/2;
  n3 = (am3 + 1) * (am3 + 2)/2;
  n4 = (am4 + 1) * (am4 + 2)/2;
  for(a=0; a<n1; a++) {
    for(b=0; b<n2; b++) {
      for(c=0; c<n3; c++) {
        for(d=0; d<n4; d++) {
          //printf("a = %d b = %d c = %d d = %d (ab|cd) = %20.15lf\n", a, b, c, d, *eri_shell_set);
          ++eri_shell_set;
        }
      }
    }
  }
}

int max(int a, int b) {
  return a >= b ? a : b;
}

void test_c_api() {
  unsigned int max_am;
  unsigned int am1, am2, am3, am4;
  double alpha1, alpha2, alpha3, alpha4;
  double A[3], B[3], C[3], D[3];

  // initialize df
  df[0] = 1;
  df[1] = 1;
  df[2] = 1;
  for(unsigned long i=3; i<maxfac*2; i++) {
    df[i] = (i-1)*df[i-2];
  }

  libint2_static_init();

  am1 = am2 = am3 = am4 = 1;
  alpha1 = 1.1;
  alpha2 = 2.3;
  alpha3 = 3.4;
  alpha4 = 4.8;
  A[0] = 0.0;  A[1] = 1.0;  A[2] = 2.0;
  B[0] = 1.0;  B[1] = 2.0;  B[2] = 0.0;
  C[0] = 2.0;  C[1] = 0.0;  C[2] = 1.0;
  D[0] = 0.0;  D[1] = 1.0;  D[2] = 2.0;

  Libint_t erieval;
  max_am = max(max(am1,am2),max(am3,am4));
  libint2_init_eri(&erieval,max_am,0);
#if LIBINT_CONTRACTED_INTS
  // if have support for contracted integrals, set the contraction length to 1
  erieval.contrdepth = 1;
#endif

  compute_eri(&erieval, am1, alpha1, A, am2, alpha2, B, am3, alpha3, C, am4, alpha4, D);

  libint2_cleanup_eri(&erieval);

  libint2_static_cleanup();
}

