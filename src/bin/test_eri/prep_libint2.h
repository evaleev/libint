#include <boost/random.hpp>

#include <cmath>
#include <vector>
#include <cstdlib>
#include <libint2.h>
#include <test_eri/eri.h>

typedef unsigned int uint;
struct RandomShellQuartetSet {
  public:
    RandomShellQuartetSet(uint* am,
                          uint veclen,
                          uint contrdepth) {

      std::copy(am, am+4, l);

      boost::mt19937 rng;                 // produces randomness out of thin air
      boost::uniform_real<> rdist(0.1, 3.0);  // distribution that maps to 0.1 .. 3.0
      boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
               die(rng, rdist);             // glues randomness with mapping


      for(uint c=0; c<4; ++c) {

        R[c].resize(3);  generate(R[c].begin(), R[c].end(), die);

        exp[c].resize(veclen);
        coef[c].resize(veclen);
        for(uint v=0; v<veclen; ++v) {
          exp[c][v].resize(contrdepth); generate(exp[c][v].begin(), exp[c][v].end(), die);
          coef[c][v].resize(contrdepth); generate(coef[c][v].begin(), coef[c][v].end(), die);
        }
      }

    }

    uint l[4];                                  // angular momenta
    std::vector<double> R[4];                   // origins
    std::vector< std::vector<double> > exp[4];  // exponents
    std::vector< std::vector<double> > coef[4]; // coefficients
};

template<typename LibintEval>
void prep_libint2(std::vector<LibintEval>& erievals,
                  const RandomShellQuartetSet& rsqset,
                  int norm_flag) {
  const uint veclen = rsqset.exp[0].size();
  const uint contrdepth = rsqset.exp[0][0].size();
  const uint contrdepth4 = erievals.size();

  const double* A = &(rsqset.R[0][0]);
  const double* B = &(rsqset.R[1][0]);
  const double* C = &(rsqset.R[2][0]);
  const double* D = &(rsqset.R[3][0]);

  const uint* am = rsqset.l;
  const unsigned int amtot = am[0] + am[1] + am[2] + am[3] + 6;
  double* F = init_array(amtot + 1);

  uint p0123 = 0;
  for (uint p0 = 0; p0 < contrdepth; p0++) {
    for (uint p1 = 0; p1 < contrdepth; p1++) {
      for (uint p2 = 0; p2 < contrdepth; p2++) {
        for (uint p3 = 0; p3 < contrdepth; p3++, p0123++) {

          LibintEval* erieval = &erievals[p0123];

          for (uint v = 0; v < veclen; v++) {

            const double alpha0 = rsqset.exp[0][v][p0];
            const double alpha1 = rsqset.exp[1][v][p1];
            const double alpha2 = rsqset.exp[2][v][p2];
            const double alpha3 = rsqset.exp[3][v][p3];

            const double c0 = rsqset.coef[0][v][p0];
            const double c1 = rsqset.coef[1][v][p1];
            const double c2 = rsqset.coef[2][v][p2];
            const double c3 = rsqset.coef[3][v][p3];

            const double gammap = alpha0 + alpha1;
            const double rhop = alpha0 * alpha1 / gammap;
            const double Px = (alpha0 * A[0] + alpha1 * B[0]) / gammap;
            const double Py = (alpha0 * A[1] + alpha1 * B[1]) / gammap;
            const double Pz = (alpha0 * A[2] + alpha1 * B[2]) / gammap;
            const double PAx = Px - A[0];
            const double PAy = Py - A[1];
            const double PAz = Pz - A[2];
            const double PBx = Px - B[0];
            const double PBy = Py - B[1];
            const double PBz = Pz - B[2];
            const double AB_x = A[0] - B[0];
            const double AB_y = A[1] - B[1];
            const double AB_z = A[2] - B[2];
            const double AB2 = AB_x * AB_x + AB_y * AB_y + AB_z * AB_z;

#if LIBINT2_DEFINED(eri,PA_x)
            erieval->PA_x[v] = PAx;
#endif
#if LIBINT2_DEFINED(eri,PA_y)
            erieval->PA_y[v] = PAy;
#endif
#if LIBINT2_DEFINED(eri,PA_z)
            erieval->PA_z[v] = PAz;
#endif
#if LIBINT2_DEFINED(eri,PB_x)
            erieval->PB_x[v] = PBx;
#endif
#if LIBINT2_DEFINED(eri,PB_y)
            erieval->PB_y[v] = PBy;
#endif
#if LIBINT2_DEFINED(eri,PB_z)
            erieval->PB_z[v] = PBz;
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

            const double gammaq = alpha2 + alpha3;
            const double rhoq = alpha2 * alpha3 / gammaq;
            const double gammapq = gammap * gammaq / (gammap + gammaq);
            const double Qx = (alpha2 * C[0] + alpha3 * D[0]) / gammaq;
            const double Qy = (alpha2 * C[1] + alpha3 * D[1]) / gammaq;
            const double Qz = (alpha2 * C[2] + alpha3 * D[2]) / gammaq;
            const double QCx = Qx - C[0];
            const double QCy = Qy - C[1];
            const double QCz = Qz - C[2];
            const double QDx = Qx - D[0];
            const double QDy = Qy - D[1];
            const double QDz = Qz - D[2];
            const double CD_x = C[0] - D[0];
            const double CD_y = C[1] - D[1];
            const double CD_z = C[2] - D[2];
            const double CD2 = CD_x * CD_x + CD_y * CD_y + CD_z * CD_z;

#if LIBINT2_DEFINED(eri,QC_x)
            erieval->QC_x[v] = QCx;
#endif
#if LIBINT2_DEFINED(eri,QC_y)
            erieval->QC_y[v] = QCy;
#endif
#if LIBINT2_DEFINED(eri,QC_z)
            erieval->QC_z[v] = QCz;
#endif
#if LIBINT2_DEFINED(eri,QD_x)
            erieval->QD_x[v] = QDx;
#endif
#if LIBINT2_DEFINED(eri,QD_y)
            erieval->QD_y[v] = QDy;
#endif
#if LIBINT2_DEFINED(eri,QD_z)
            erieval->QD_z[v] = QDz;
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

    // Prefactors for interelectron transfer relation
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_x)
            erieval->TwoPRepITR_pfac0_0_x[v] = - (alpha1*AB_x + alpha3*CD_x)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_y)
            erieval->TwoPRepITR_pfac0_0_y[v] = - (alpha1*AB_y + alpha3*CD_y)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_z)
            erieval->TwoPRepITR_pfac0_0_z[v] = - (alpha1*AB_z + alpha3*CD_z)/gammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_x)
            erieval->TwoPRepITR_pfac0_1_x[v] = - (alpha1*AB_x + alpha3*CD_x)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_y)
            erieval->TwoPRepITR_pfac0_1_y[v] = - (alpha1*AB_y + alpha3*CD_y)/gammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_z)
            erieval->TwoPRepITR_pfac0_1_z[v] = - (alpha1*AB_z + alpha3*CD_z)/gammaq;
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
            const double PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;
            const double Wx = (gammap * Px + gammaq * Qx) / (gammap + gammaq);
            const double Wy = (gammap * Py + gammaq * Qy) / (gammap + gammaq);
            const double Wz = (gammap * Pz + gammaq * Qz) / (gammap + gammaq);

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

            // prefactors for derivative ERI relations
#if LIBINT2_DEFINED(eri,alpha1_rho_over_zeta2)
            erieval->alpha1_rho_over_zeta2[v] = alpha0 * gammapq / (gammap * gammap);
#endif
#if LIBINT2_DEFINED(eri,alpha2_rho_over_zeta2)
            erieval->alpha2_rho_over_zeta2[v] = alpha1 * gammapq / (gammap * gammap);
#endif
#if LIBINT2_DEFINED(eri,alpha3_rho_over_eta2)
            erieval->alpha3_rho_over_eta2[v] = alpha2 * gammapq / (gammaq * gammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha4_rho_over_eta2)
            erieval->alpha4_rho_over_eta2[v] = alpha3 * gammapq / (gammaq * gammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha1_over_zetapluseta)
            erieval->alpha1_over_zetapluseta[v] = alpha0 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha2_over_zetapluseta)
            erieval->alpha2_over_zetapluseta[v] = alpha1 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha3_over_zetapluseta)
            erieval->alpha3_over_zetapluseta[v] = alpha2 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha4_over_zetapluseta)
            erieval->alpha4_over_zetapluseta[v] = alpha3 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri,rho12_over_alpha1)
            erieval->rho12_over_alpha1[v] = rhop / alpha0;
#endif
#if LIBINT2_DEFINED(eri,rho12_over_alpha2)
            erieval->rho12_over_alpha2[v] = rhop / alpha1;
#endif
#if LIBINT2_DEFINED(eri,rho34_over_alpha3)
            erieval->rho34_over_alpha3[v] = rhoq / alpha2;
#endif
#if LIBINT2_DEFINED(eri,rho34_over_alpha4)
            erieval->rho34_over_alpha4[v] = rhoq / alpha3;
#endif
#if LIBINT2_DEFINED(eri,two_alpha0_bra)
            erieval->two_alpha0_bra[v] = 2.0 * alpha0;
#endif
#if LIBINT2_DEFINED(eri,two_alpha0_ket)
            erieval->two_alpha0_ket[v] = 2.0 * alpha1;
#endif
#if LIBINT2_DEFINED(eri,two_alpha1_bra)
            erieval->two_alpha1_bra[v] = 2.0 * alpha2;
#endif
#if LIBINT2_DEFINED(eri,two_alpha1_ket)
            erieval->two_alpha1_ket[v] = 2.0 * alpha3;
#endif

            double K1 = exp(-alpha0 * alpha1 * AB2 / gammap);
            double K2 = exp(-alpha2 * alpha3 * CD2 / gammaq);
            double pfac = 2 * pow(M_PI, 2.5) * K1 * K2 / (gammap * gammaq * sqrt(gammap
                                                                                 + gammaq));
            pfac *= c0 * c1 * c2 * c3;

            if (norm_flag > 0) {
//              pfac *= norm_const(l1,m1,n1,alpha0,A);
//              pfac *= norm_const(l2,m2,n2,alpha1,B);
//              pfac *= norm_const(l3,m3,n3,alpha2,C);
//              pfac *= norm_const(l4,m4,n4,alpha3,D);
            }

            calc_f(F, amtot, PQ2 * gammapq);

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
      }
    }
  } // end of primitive loops

  free_array(F);
}

