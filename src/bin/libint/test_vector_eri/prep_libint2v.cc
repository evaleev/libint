
#include <cmath>
#include <stdlib.h>
#include <libint2.h>
#include <eri.h>

#define MAXAM 8
#define MAXFAC 100
#define EPS 1.0E-17     /* Absolute precision in computing Fm(t)
                           (see recursion:calc_fij() ) */
#define MIN(a,b) ((a)>(b) ? (b) : (a))

void prep_libint2v(Libint_t* libint, unsigned int veclength, unsigned int am1,
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

  const double PQx = Px - Qx;
  const double PQy = Py - Qy;
  const double PQz = Pz - Qz;
  const double PQ2 = PQx*PQx + PQy*PQy + PQz*PQz;
  const double Wx = (gammap*Px + gammaq*Qx)/(gammap+gammaq);
  const double Wy = (gammap*Py + gammaq*Qy)/(gammap+gammaq);
  const double Wz = (gammap*Pz + gammaq*Qz)/(gammap+gammaq);

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


  libint->veclength = veclength;
  libint->PA_x = new double[veclength];
  libint->PA_y = new double[veclength];
  libint->PA_z = new double[veclength];
  libint->QC_x = new double[veclength];
  libint->QC_y = new double[veclength];
  libint->QC_z = new double[veclength];
  libint->AB_x = new double[veclength];
  libint->AB_y = new double[veclength];
  libint->AB_z = new double[veclength];
  libint->CD_x = new double[veclength];
  libint->CD_y = new double[veclength];
  libint->CD_z = new double[veclength];
  libint->WP_x = new double[veclength];
  libint->WP_y = new double[veclength];
  libint->WP_z = new double[veclength];
  libint->WQ_x = new double[veclength];
  libint->WQ_y = new double[veclength];
  libint->WQ_z = new double[veclength];
  libint->oo2z = new double[veclength];
  libint->oo2e = new double[veclength];
  libint->oo2ze = new double[veclength];
  libint->roz = new double[veclength];
  libint->roe = new double[veclength];
  libint->__ss_1_over_r_12_ss___up_0 = new double[veclength];
  libint->__ss_1_over_r_12_ss___up_1 = new double[veclength];
  libint->__ss_1_over_r_12_ss___up_2 = new double[veclength];
  libint->__ss_1_over_r_12_ss___up_3 = new double[veclength];
  libint->__ss_1_over_r_12_ss___up_4 = new double[veclength];
  libint->__ss_1_over_r_12_ss___up_5 = new double[veclength];
  libint->__ss_1_over_r_12_ss___up_6 = new double[veclength];
  libint->__ss_1_over_r_12_ss___up_7 = new double[veclength];
  libint->__ss_1_over_r_12_ss___up_8 = new double[veclength];
  libint->__ss_1_over_r_12_ss___up_9 = new double[veclength];
  libint->__ss_1_over_r_12_ss___up_10 = new double[veclength];
  libint->__ss_1_over_r_12_ss___up_11 = new double[veclength];
  libint->__ss_1_over_r_12_ss___up_12 = new double[veclength];

  for(int v=0; v<veclength; v++) {
  
    libint->PA_x[v] = PAx;
    libint->PA_y[v] = PAy;
    libint->PA_z[v] = PAz;
    libint->AB_x[v] = A[0] - B[0];
    libint->AB_y[v] = A[1] - B[1];
    libint->AB_z[v] = A[2] - B[2];
    libint->oo2z[v] = 0.5/gammap;
    
    libint->QC_x[v] = QCx;
    libint->QC_y[v] = QCy;
    libint->QC_z[v] = QCz;
    libint->CD_x[v] = C[0] - D[0];
    libint->CD_y[v] = C[1] - D[1];
    libint->CD_z[v] = C[2] - D[2];
    libint->oo2e[v] = 0.5/gammaq;
    
    libint->WP_x[v] = Wx - Px;
    libint->WP_y[v] = Wy - Py;
    libint->WP_z[v] = Wz - Pz;
    libint->WQ_x[v] = Wx - Qx;
    libint->WQ_y[v] = Wy - Qy;
    libint->WQ_z[v] = Wz - Qz;
    libint->oo2ze[v] = 0.5/(gammap+gammaq);
    libint->roz[v] = gammapq/gammap;
    libint->roe[v] = gammapq/gammaq;
  
    libint->__ss_1_over_r_12_ss___up_0[v] = pfac*F[0];
    libint->__ss_1_over_r_12_ss___up_1[v] = pfac*F[1];
    libint->__ss_1_over_r_12_ss___up_2[v] = pfac*F[2];
    libint->__ss_1_over_r_12_ss___up_3[v] = pfac*F[3];
    libint->__ss_1_over_r_12_ss___up_4[v] = pfac*F[4];
    libint->__ss_1_over_r_12_ss___up_5[v] = pfac*F[5];
    libint->__ss_1_over_r_12_ss___up_6[v] = pfac*F[6];
    libint->__ss_1_over_r_12_ss___up_7[v] = pfac*F[7];
    libint->__ss_1_over_r_12_ss___up_8[v] = pfac*F[8];
    libint->__ss_1_over_r_12_ss___up_9[v] = pfac*F[9];
    libint->__ss_1_over_r_12_ss___up_10[v] = pfac*F[10];
    libint->__ss_1_over_r_12_ss___up_11[v] = pfac*F[11];
    libint->__ss_1_over_r_12_ss___up_12[v] = pfac*F[12];
  }
  
}



