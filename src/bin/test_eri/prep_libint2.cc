
#include <cmath>
#include <vector>
#include <cstdlib>
#include <libint2.h>
#include <test_eri/eri.h>


#define MAXAM 8
#define MAXFAC 100
#define EPS 1.0E-17     /* Absolute precision in computing Fm(t)
                           (see recursion:calc_fij() ) */
#define MIN(a,b) ((a)>(b) ? (b) : (a))

void prep_libint2(Libint_t* libint, unsigned int am1,
		  const std::vector<double>& alpha1, double A[3],
		  unsigned int am2, 
		  const std::vector<double>& alpha2, double B[3],
		  unsigned int am3, 
		  const std::vector<double>& alpha3, double C[3],
		  unsigned int am4, 
		  const std::vector<double>& alpha4, double D[3], int norm_flag,
                  unsigned int veclen)
{
  for(int v=0; v<veclen; v++) {
    const double gammap = alpha1[v] + alpha2[v];
    const double Px = (alpha1[v]*A[0] + alpha2[v]*B[0])/gammap;
    const double Py = (alpha1[v]*A[1] + alpha2[v]*B[1])/gammap;
    const double Pz = (alpha1[v]*A[2] + alpha2[v]*B[2])/gammap;
    const double PAx = Px - A[0];
    const double PAy = Py - A[1];
    const double PAz = Pz - A[2];
    const double PBx = Px - B[0];
    const double PBy = Py - B[1];
    const double PBz = Pz - B[2];
    const double AB2 = (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) + (A[2]-B[2])*(A[2]-B[2]);
    
    libint->PA_x[v] = PAx;
    libint->PA_y[v] = PAy;
    libint->PA_z[v] = PAz;
    libint->AB_x[v] = A[0] - B[0];
    libint->AB_y[v] = A[1] - B[1];
    libint->AB_z[v] = A[2] - B[2];
    libint->oo2z[v] = 0.5/gammap;
    
    const double gammaq = alpha3[v] + alpha4[v];
    const double gammapq = gammap*gammaq/(gammap+gammaq);
    const double Qx = (alpha3[v]*C[0] + alpha4[v]*D[0])/gammaq;
    const double Qy = (alpha3[v]*C[1] + alpha4[v]*D[1])/gammaq;
    const double Qz = (alpha3[v]*C[2] + alpha4[v]*D[2])/gammaq;
    const double QCx = Qx - C[0];
    const double QCy = Qy - C[1];
    const double QCz = Qz - C[2];
    const double QDx = Qx - D[0];
    const double QDy = Qy - D[1];
    const double QDz = Qz - D[2];
    const double CD2 = (C[0]-D[0])*(C[0]-D[0]) + (C[1]-D[1])*(C[1]-D[1]) + (C[2]-D[2])*(C[2]-D[2]);
    
    libint->QC_x[v] = QCx;
    libint->QC_y[v] = QCy;
    libint->QC_z[v] = QCz;
    libint->CD_x[v] = C[0] - D[0];
    libint->CD_y[v] = C[1] - D[1];
    libint->CD_z[v] = C[2] - D[2];
    libint->oo2e[v] = 0.5/gammaq;
    
    const double PQx = Px - Qx;
    const double PQy = Py - Qy;
    const double PQz = Pz - Qz;
    const double PQ2 = PQx*PQx + PQy*PQy + PQz*PQz;
    const double Wx = (gammap*Px + gammaq*Qx)/(gammap+gammaq);
    const double Wy = (gammap*Py + gammaq*Qy)/(gammap+gammaq);
    const double Wz = (gammap*Pz + gammaq*Qz)/(gammap+gammaq);
    
    libint->WP_x[v] = Wx - Px;
    libint->WP_y[v] = Wy - Py;
    libint->WP_z[v] = Wz - Pz;
    libint->WQ_x[v] = Wx - Qx;
    libint->WQ_y[v] = Wy - Qy;
    libint->WQ_z[v] = Wz - Qz;
    libint->oo2ze[v] = 0.5/(gammap+gammaq);
    libint->roz[v] = gammapq/gammap;
    libint->roe[v] = gammapq/gammaq;
    
    double K1 = exp(-alpha1[v]*alpha2[v]*AB2/gammap);
    double K2 = exp(-alpha3[v]*alpha4[v]*CD2/gammaq);
    double pfac = 2*pow(M_PI,2.5)*K1*K2/(gammap*gammaq*sqrt(gammap+gammaq));
    
    if (norm_flag > 0) {
      /*    pfac *= norm_const(l1,m1,n1,alpha1[v],A);
      pfac *= norm_const(l2,m2,n2,alpha2[v],B);
      pfac *= norm_const(l3,m3,n3,alpha3[v],C);
      pfac *= norm_const(l4,m4,n4,alpha4[v],D);*/
    }
    
    unsigned int am = am1 + am2 + am3 + am4;
    double* F = init_array(am+1);
    calc_f(F,am,PQ2*gammapq);
    
    // using dangerous macros from libint2.h
    libint->LIBINT_T_SS_EREP_SS(0)[v] = pfac*F[0];
    libint->LIBINT_T_SS_EREP_SS(1)[v] = pfac*F[1];
    libint->LIBINT_T_SS_EREP_SS(2)[v] = pfac*F[2];
    libint->LIBINT_T_SS_EREP_SS(3)[v] = pfac*F[3];
    libint->LIBINT_T_SS_EREP_SS(4)[v] = pfac*F[4];
    libint->LIBINT_T_SS_EREP_SS(5)[v] = pfac*F[5];
    libint->LIBINT_T_SS_EREP_SS(6)[v] = pfac*F[6];
    libint->LIBINT_T_SS_EREP_SS(7)[v] = pfac*F[7];
    libint->LIBINT_T_SS_EREP_SS(8)[v] = pfac*F[8];
    libint->LIBINT_T_SS_EREP_SS(9)[v] = pfac*F[9];
    libint->LIBINT_T_SS_EREP_SS(10)[v] = pfac*F[10];
    libint->LIBINT_T_SS_EREP_SS(11)[v] = pfac*F[11];
    libint->LIBINT_T_SS_EREP_SS(12)[v] = pfac*F[12];
    libint->LIBINT_T_SS_EREP_SS(13)[v] = pfac*F[13];
    libint->LIBINT_T_SS_EREP_SS(14)[v] = pfac*F[14];
    libint->LIBINT_T_SS_EREP_SS(15)[v] = pfac*F[15];
    libint->LIBINT_T_SS_EREP_SS(16)[v] = pfac*F[16];
    libint->LIBINT_T_SS_EREP_SS(17)[v] = pfac*F[17];
    libint->LIBINT_T_SS_EREP_SS(18)[v] = pfac*F[18];
    libint->LIBINT_T_SS_EREP_SS(19)[v] = pfac*F[19];
    libint->LIBINT_T_SS_EREP_SS(20)[v] = pfac*F[20];
  } // end of v loop
}



