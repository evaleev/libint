
#include <cmath>
#include <stdlib.h>
#include <libint2.h>
#include <test_eri/eri.h>


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
  
  libint->PA_x[0] = PAx;
  libint->PA_y[0] = PAy;
  libint->PA_z[0] = PAz;
  libint->AB_x[0] = A[0] - B[0];
  libint->AB_y[0] = A[1] - B[1];
  libint->AB_z[0] = A[2] - B[2];
  libint->oo2z[0] = 0.5/gammap;
  
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

  libint->QC_x[0] = QCx;
  libint->QC_y[0] = QCy;
  libint->QC_z[0] = QCz;
  libint->CD_x[0] = C[0] - D[0];
  libint->CD_y[0] = C[1] - D[1];
  libint->CD_z[0] = C[2] - D[2];
  libint->oo2e[0] = 0.5/gammaq;

  const double PQx = Px - Qx;
  const double PQy = Py - Qy;
  const double PQz = Pz - Qz;
  const double PQ2 = PQx*PQx + PQy*PQy + PQz*PQz;
  const double Wx = (gammap*Px + gammaq*Qx)/(gammap+gammaq);
  const double Wy = (gammap*Py + gammaq*Qy)/(gammap+gammaq);
  const double Wz = (gammap*Pz + gammaq*Qz)/(gammap+gammaq);

  libint->WP_x[0] = Wx - Px;
  libint->WP_y[0] = Wy - Py;
  libint->WP_z[0] = Wz - Pz;
  libint->WQ_x[0] = Wx - Qx;
  libint->WQ_y[0] = Wy - Qy;
  libint->WQ_z[0] = Wz - Qz;
  libint->oo2ze[0] = 0.5/(gammap+gammaq);
  libint->roz[0] = gammapq/gammap;
  libint->roe[0] = gammapq/gammaq;
  
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

  // using dangerous macros from libint2.h
  libint->LIBINT_T_SS_EREP_SS(0)[0] = pfac*F[0];
  libint->LIBINT_T_SS_EREP_SS(1)[0] = pfac*F[1];
  libint->LIBINT_T_SS_EREP_SS(2)[0] = pfac*F[2];
  libint->LIBINT_T_SS_EREP_SS(3)[0] = pfac*F[3];
  libint->LIBINT_T_SS_EREP_SS(4)[0] = pfac*F[4];
  libint->LIBINT_T_SS_EREP_SS(5)[0] = pfac*F[5];
  libint->LIBINT_T_SS_EREP_SS(6)[0] = pfac*F[6];
  libint->LIBINT_T_SS_EREP_SS(7)[0] = pfac*F[7];
  libint->LIBINT_T_SS_EREP_SS(8)[0] = pfac*F[8];
  libint->LIBINT_T_SS_EREP_SS(9)[0] = pfac*F[9];
  libint->LIBINT_T_SS_EREP_SS(10)[0] = pfac*F[10];
  libint->LIBINT_T_SS_EREP_SS(11)[0] = pfac*F[11];
  libint->LIBINT_T_SS_EREP_SS(12)[0] = pfac*F[12];
  libint->LIBINT_T_SS_EREP_SS(13)[0] = pfac*F[13];
  libint->LIBINT_T_SS_EREP_SS(14)[0] = pfac*F[14];
  libint->LIBINT_T_SS_EREP_SS(15)[0] = pfac*F[15];
  libint->LIBINT_T_SS_EREP_SS(16)[0] = pfac*F[16];
  libint->LIBINT_T_SS_EREP_SS(17)[0] = pfac*F[17];
  libint->LIBINT_T_SS_EREP_SS(18)[0] = pfac*F[18];
  libint->LIBINT_T_SS_EREP_SS(19)[0] = pfac*F[19];
  libint->LIBINT_T_SS_EREP_SS(20)[0] = pfac*F[20];
  
}



