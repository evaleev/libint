#include <iostream>
#include <algorithm>

#include <libint2.h>
#include <libint2/boys.h>

using namespace std;
using namespace libint2;

/** This function evaluates ERI over 4 primitive Gaussian shells.
    See tests/eri/test.cc for an example of how to deal with
    contracted Gaussians.

    For simplicity, many details are omitted here, e.g. normalization.
  */
    void
    compute_eri(unsigned int am1, double alpha1, double A[3],
                unsigned int am2, double alpha2, double B[3],
                unsigned int am3, double alpha3, double C[3],
                unsigned int am4, double alpha4, double D[3]
               )
{
    // I will assume that libint2_static_init() has been called elsewhere!

    //-------------------------------------------------------------------------------------
    // ------ the code in this section would usually be placed outside this function ------
    //        to occur once per calculation, not every time this function is called
    // - allocate ERI evaluator object
    Libint_eri_t erieval;
    const unsigned int max_am = max(max(am1,am2),max(am3,am4));
    libint2_init_eri(&erieval,max_am,0);
#if LIBINT_CONTRACTED_INTS
    // if have support for contracted integrals, set the contraction length to 1
    erieval.contrdepth = 1;
#endif
    // - initialize Boys function evaluator
    FmEval_Chebyshev7 fmeval(max_am);
    //-------------------------------------------------------------------------------------

    //
    // Compute requisite data -- many of these quantities would be precomputed
    // for all nonnegligible shell pairs somewhere else
    //
    const double gammap = alpha1 + alpha2;
    const double Px = (alpha1*A[0] + alpha2*B[0])/gammap;
    const double Py = (alpha1*A[1] + alpha2*B[1])/gammap;
    const double Pz = (alpha1*A[2] + alpha2*B[2])/gammap;
    const double PAx = Px - A[0];
    const double PAy = Py - A[1];
    const double PAz = Pz - A[2];
    const double PBx = Px - B[0];
    const double PBy = Py - B[1];
    const double PBz = Pz - B[2];
    const double AB2 = (A[0]-B[0])*(A[0]-B[0])
                     + (A[1]-B[1])*(A[1]-B[1])
                     + (A[2]-B[2])*(A[2]-B[2]);

    erieval.PA_x[0] = PAx;
    erieval.PA_y[0] = PAy;
    erieval.PA_z[0] = PAz;
    erieval.AB_x[0] = A[0] - B[0];
    erieval.AB_y[0] = A[1] - B[1];
    erieval.AB_z[0] = A[2] - B[2];
    erieval.oo2z[0] = 0.5/gammap;

    const double gammaq = alpha3 + alpha4;
    const double gammapq = gammap*gammaq/(gammap+gammaq);
    const double Qx = (alpha3*C[0] + alpha4*D[0])/gammaq;
    const double Qy = (alpha3*C[1] + alpha4*D[1])/gammaq;
    const double Qz = (alpha3*C[2] + alpha4*D[2])/gammaq;
    const double QCx = Qx - C[0];
    const double QCy = Qy - C[1];
    const double QCz = Qz - C[2];
    const double QDx = Qx - D[0];
    const double QDy = Qy - D[1];
    const double QDz = Qz - D[2];
    const double CD2 = (C[0]-D[0])*(C[0]-D[0])
                     + (C[1]-D[1])*(C[1]-D[1])
                     + (C[2]-D[2])*(C[2]-D[2]);

    erieval.QC_x[0] = QCx;
    erieval.QC_y[0] = QCy;
    erieval.QC_z[0] = QCz;
    erieval.CD_x[0] = C[0] - D[0];
    erieval.CD_y[0] = C[1] - D[1];
    erieval.CD_z[0] = C[2] - D[2];
    erieval.oo2e[0] = 0.5/gammaq;

    const double PQx = Px - Qx;
    const double PQy = Py - Qy;
    const double PQz = Pz - Qz;
    const double PQ2 = PQx*PQx + PQy*PQy + PQz*PQz;
    const double Wx = (gammap*Px + gammaq*Qx)/(gammap+gammaq);
    const double Wy = (gammap*Py + gammaq*Qy)/(gammap+gammaq);
    const double Wz = (gammap*Pz + gammaq*Qz)/(gammap+gammaq);

    erieval.WP_x[0] = Wx - Px;
    erieval.WP_y[0] = Wy - Py;
    erieval.WP_z[0] = Wz - Pz;
    erieval.WQ_x[0] = Wx - Qx;
    erieval.WQ_y[0] = Wy - Qy;
    erieval.WQ_z[0] = Wz - Qz;
    erieval.oo2ze[0] = 0.5/(gammap+gammaq);
    erieval.roz[0] = gammapq/gammap;
    erieval.roe[0] = gammapq/gammaq;

    double K1 = exp(-alpha1*alpha2*AB2/gammap);
    double K2 = exp(-alpha3*alpha4*CD2/gammaq);
    double pfac = 2*pow(M_PI,2.5)*K1*K2/(gammap*gammaq*sqrt(gammap+gammaq));

    //
    // evaluate Boys function F_m for all m in [0,am]
    //
    unsigned int am = am1 + am2 + am3 + am4;
    double* F = new double[am+1];
    fmeval.eval(F, PQ2*gammapq, am);

    // (00|00)^m = pfac * F_m
    erieval.LIBINT_T_SS_EREP_SS(0)[0] = pfac*F[0];
    erieval.LIBINT_T_SS_EREP_SS(1)[0] = pfac*F[1];
    erieval.LIBINT_T_SS_EREP_SS(2)[0] = pfac*F[2];
    erieval.LIBINT_T_SS_EREP_SS(3)[0] = pfac*F[3];
    erieval.LIBINT_T_SS_EREP_SS(4)[0] = pfac*F[4];
    // etc.

    // compute ERIs
    libint2_build_eri[am1][am2][am3][am4](&erieval);

    // Print out the integrals
    const double* eri_shell_set = erieval.targets[0];
    const unsigned int n1 = (am1 + 1) * (am1 + 2)/2;
    const unsigned int n2 = (am2 + 1) * (am2 + 2)/2;
    const unsigned int n3 = (am3 + 1) * (am3 + 2)/2;
    const unsigned int n4 = (am4 + 1) * (am4 + 2)/2;
    for(int a=0; a<n1; a++) {
      for(int b=0; b<n2; b++) {
        for(int c=0; c<n3; c++) {
          for(int d=0; d<n4; d++) {
            cout << "a = " << a
                 << "b = " << b
                 << "c = " << c
                 << "d = " << d
                 << "(ab|cd) = " << *eri_shell_set;
            ++eri_shell_set;
          }
        }
      }
    }

    // ------ like the code at the beginning, this usually goes outside this function ------
    libint2_cleanup_eri(&erieval);
}
