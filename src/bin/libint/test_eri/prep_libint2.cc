
#include <cmath>
#include <stdlib.h>
#include <libint2.h>
#include <eri.h>


#define MAXAM 8
#define MAXFAC 100
#define EPS 1.0E-17     /* Absolute precision in computing Fm(t)
                           (see recursion:calc_fij() ) */
#define MIN(a,b) ((a)>(b) ? (b) : (a))

void prep_libint2(Libint_t* libint, unsigned int am1,
           double alpha1, double A[3],
	   unsigned int am2, 
           double alpha2, double B[3],
	   unsigned int am3, 
           double alpha3, double C[3],
	   unsigned int am4, 
           double alpha4, double D[3], int norm_flag)
{

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
  const double AB2 = (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) + (A[2]-B[2])*(A[2]-B[2]);
  
  libint->PA_x = PAx;
  libint->PA_y = PAy;
  libint->PA_z = PAz;
  libint->AB_x = A[0] - B[0];
  libint->AB_y = A[1] - B[1];
  libint->AB_z = A[2] - B[2];
  libint->oo2z = 0.5/gammap;
  
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
  const double CD2 = (C[0]-D[0])*(C[0]-D[0]) + (C[1]-D[1])*(C[1]-D[1]) + (C[2]-D[2])*(C[2]-D[2]);

  libint->QC_x = QCx;
  libint->QC_y = QCy;
  libint->QC_z = QCz;
  libint->CD_x = C[0] - D[0];
  libint->CD_y = C[1] - D[1];
  libint->CD_z = C[2] - D[2];
  libint->oo2e = 0.5/gammaq;

  const double PQx = Px - Qx;
  const double PQy = Py - Qy;
  const double PQz = Pz - Qz;
  const double PQ2 = PQx*PQx + PQy*PQy + PQz*PQz;
  const double Wx = (gammap*Px + gammaq*Qx)/(gammap+gammaq);
  const double Wy = (gammap*Py + gammaq*Qy)/(gammap+gammaq);
  const double Wz = (gammap*Pz + gammaq*Qz)/(gammap+gammaq);

  libint->WP_x = Wx - Px;
  libint->WP_y = Wy - Py;
  libint->WP_z = Wz - Pz;
  libint->WQ_x = Wx - Qx;
  libint->WQ_y = Wy - Qy;
  libint->WQ_z = Wz - Qz;
  libint->oo2ze = 0.5/(gammap+gammaq);
  libint->roz = gammapq/gammap;
  libint->roe = gammapq/gammaq;
  
  double K1 = exp(-alpha1*alpha2*AB2/gammap);
  double K2 = exp(-alpha3*alpha4*CD2/gammaq);
  double pfac = 2*pow(M_PI,2.5)*K1*K2/(gammap*gammaq*sqrt(gammap+gammaq));

  if (norm_flag > 0) {
/*    pfac *= norm_const(l1,m1,n1,alpha1,A);
    pfac *= norm_const(l2,m2,n2,alpha2,B);
    pfac *= norm_const(l3,m3,n3,alpha3,C);
    pfac *= norm_const(l4,m4,n4,alpha4,D);*/
  }
  
  unsigned int am = am1 + am2 + am3 + am4;
  double* F = init_array(am+1);
  calc_f(F,am,PQ2*gammapq);
  
  libint->__ss_1_over_r_12_ss___up_0 = pfac*F[0];
  libint->__ss_1_over_r_12_ss___up_1 = pfac*F[1];
  libint->__ss_1_over_r_12_ss___up_2 = pfac*F[2];
  libint->__ss_1_over_r_12_ss___up_3 = pfac*F[3];
  libint->__ss_1_over_r_12_ss___up_4 = pfac*F[4];
  libint->__ss_1_over_r_12_ss___up_5 = pfac*F[5];
  libint->__ss_1_over_r_12_ss___up_6 = pfac*F[6];
  libint->__ss_1_over_r_12_ss___up_7 = pfac*F[7];
  libint->__ss_1_over_r_12_ss___up_8 = pfac*F[8];
  libint->__ss_1_over_r_12_ss___up_9 = pfac*F[9];
  libint->__ss_1_over_r_12_ss___up_10 = pfac*F[10];
  libint->__ss_1_over_r_12_ss___up_11 = pfac*F[11];
  libint->__ss_1_over_r_12_ss___up_12 = pfac*F[12];
  
}



