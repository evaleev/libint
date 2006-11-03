
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

void prep_libint2(Libint_eri_t* erieval, unsigned int am1,
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
    const double AB_x = A[0] - B[0];
    const double AB_y = A[1] - B[1];
    const double AB_z = A[2] - B[2];
    const double AB2 = AB_x*AB_x + AB_y*AB_y + AB_z*AB_z;

#define LIBINT2_DEFINED(taskname,symbol) __prescanned_libint2_defined__(taskname,symbol)
#define __prescanned_libint2_defined__(taskname,symbol) LIBINT2_DEFINED_##symbol##_##taskname

#if LIBINT2_DEFINED(eri,PA_x)
    erieval->PA_x[v] = PAx;
#endif
#if LIBINT2_DEFINED(eri,PA_y)
    erieval->PA_y[v] = PAy;
#endif
#if LIBINT2_DEFINED(eri,PA_z)
    erieval->PA_z[v] = PAz;
#endif

#if LIBINT2_DEFINED(eri,AB_x)
    erieval->AB_x[v] = AB_x;
#endif
#if LIBINT2_DEFINED(eri,AB_y)
    erieval->AB_y[v] = AB_y;
#endif
#if LIBINT2_DEFINED(eri,AB_z)
    erieval->AB_z[v] = AB_z;
#endif
#if LIBINT2_DEFINED(eri,oo2z)
    erieval->oo2z[v] = 0.5/gammap;
#endif
    
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
    const double CD_x = C[0] - D[0];
    const double CD_y = C[1] - D[1];
    const double CD_z = C[2] - D[2];
    const double CD2 = CD_x*CD_x + CD_y*CD_y + CD_z*CD_z;
    
#if LIBINT2_DEFINED(eri,QC_x)
    erieval->QC_x[v] = QCx;
#endif
#if LIBINT2_DEFINED(eri,QC_y)
    erieval->QC_y[v] = QCy;
#endif
#if LIBINT2_DEFINED(eri,QC_z)
    erieval->QC_z[v] = QCz;
#endif

#if LIBINT2_DEFINED(eri,CD_x)
    erieval->CD_x[v] = CD_x;
#endif
#if LIBINT2_DEFINED(eri,CD_y)
    erieval->CD_y[v] = CD_y;
#endif 
#if LIBINT2_DEFINED(eri,CD_z)
    erieval->CD_z[v] = CD_z;
#endif
#if LIBINT2_DEFINED(eri,oo2e)
    erieval->oo2e[v] = 0.5/gammaq;
#endif
    
    // Prefactors for interelecttron transfer relation
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_x)
    erieval->TwoPRepITR_pfac0_0_x[v] = - (alpha2[v]*AB_x + alpha4[v]*CD_x)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_y)
    erieval->TwoPRepITR_pfac0_0_y[v] = - (alpha2[v]*AB_y + alpha4[v]*CD_y)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_z)
    erieval->TwoPRepITR_pfac0_0_z[v] = - (alpha2[v]*AB_z + alpha4[v]*CD_z)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_x)
    erieval->TwoPRepITR_pfac0_1_x[v] = - (alpha2[v]*AB_x + alpha4[v]*CD_x)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_y)
    erieval->TwoPRepITR_pfac0_1_y[v] = - (alpha2[v]*AB_y + alpha4[v]*CD_y)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_z)
    erieval->TwoPRepITR_pfac0_1_z[v] = - (alpha2[v]*AB_z + alpha4[v]*CD_z)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac1_0)
    erieval->TwoPRepITR_pfac1_0[v] = -gammaq/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac1_1)
    erieval->TwoPRepITR_pfac1_1[v] = -gammap/gammaq;
#endif
    
    const double PQx = Px - Qx;
    const double PQy = Py - Qy;
    const double PQz = Pz - Qz;
    const double PQ2 = PQx*PQx + PQy*PQy + PQz*PQz;
    const double Wx = (gammap*Px + gammaq*Qx)/(gammap+gammaq);
    const double Wy = (gammap*Py + gammaq*Qy)/(gammap+gammaq);
    const double Wz = (gammap*Pz + gammaq*Qz)/(gammap+gammaq);
    
#if LIBINT2_DEFINED(eri,WP_x)
    erieval->WP_x[v] = Wx - Px;
#endif
#if LIBINT2_DEFINED(eri,WP_y)
    erieval->WP_y[v] = Wy - Py;
#endif
#if LIBINT2_DEFINED(eri,WP_z)
    erieval->WP_z[v] = Wz - Pz;
#endif
#if LIBINT2_DEFINED(eri,WQ_x)
    erieval->WQ_x[v] = Wx - Qx;
#endif
#if LIBINT2_DEFINED(eri,WQ_y)
    erieval->WQ_y[v] = Wy - Qy;
#endif
#if LIBINT2_DEFINED(eri,WQ_z)
    erieval->WQ_z[v] = Wz - Qz;
#endif
#if LIBINT2_DEFINED(eri,oo2ze)
    erieval->oo2ze[v] = 0.5/(gammap+gammaq);
#endif
#if LIBINT2_DEFINED(eri,roz)
    erieval->roz[v] = gammapq/gammap;
#endif
#if LIBINT2_DEFINED(eri,roe)
    erieval->roe[v] = gammapq/gammaq;
#endif
    
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
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(0))
    erieval->LIBINT_T_SS_EREP_SS(0)[v] = pfac*F[0];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(1))
    erieval->LIBINT_T_SS_EREP_SS(1)[v] = pfac*F[1];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(2))
    erieval->LIBINT_T_SS_EREP_SS(2)[v] = pfac*F[2];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(3))
    erieval->LIBINT_T_SS_EREP_SS(3)[v] = pfac*F[3];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(4))
    erieval->LIBINT_T_SS_EREP_SS(4)[v] = pfac*F[4];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(5))
    erieval->LIBINT_T_SS_EREP_SS(5)[v] = pfac*F[5];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(6))
    erieval->LIBINT_T_SS_EREP_SS(6)[v] = pfac*F[6];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(7))
    erieval->LIBINT_T_SS_EREP_SS(7)[v] = pfac*F[7];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(8))
    erieval->LIBINT_T_SS_EREP_SS(8)[v] = pfac*F[8];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(9))
    erieval->LIBINT_T_SS_EREP_SS(9)[v] = pfac*F[9];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(10))
    erieval->LIBINT_T_SS_EREP_SS(10)[v] = pfac*F[10];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(11))
    erieval->LIBINT_T_SS_EREP_SS(11)[v] = pfac*F[11];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(12))
    erieval->LIBINT_T_SS_EREP_SS(12)[v] = pfac*F[12];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(13))
    erieval->LIBINT_T_SS_EREP_SS(13)[v] = pfac*F[13];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(14))
    erieval->LIBINT_T_SS_EREP_SS(14)[v] = pfac*F[14];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(15))
    erieval->LIBINT_T_SS_EREP_SS(15)[v] = pfac*F[15];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(16))
    erieval->LIBINT_T_SS_EREP_SS(16)[v] = pfac*F[16];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(17))
    erieval->LIBINT_T_SS_EREP_SS(17)[v] = pfac*F[17];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(18))
    erieval->LIBINT_T_SS_EREP_SS(18)[v] = pfac*F[18];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(19))
    erieval->LIBINT_T_SS_EREP_SS(19)[v] = pfac*F[19];
#endif
#if LIBINT2_DEFINED(eri,LIBINT_T_SS_EREP_SS(20))
    erieval->LIBINT_T_SS_EREP_SS(20)[v] = pfac*F[20];
#endif
  } // end of v loop
}



