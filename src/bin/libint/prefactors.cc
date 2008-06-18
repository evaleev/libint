
#include <smart_ptr.h>
#include <prefactors.h>

using namespace libint2;

Prefactors::Prefactors() :
  rho(new rdouble("rho")),
  one_o_2alphasum(new rdouble("oo2ze"))
{
  typedef SafePtr<rdouble> rdptr;
  typedef SafePtr<rdouble> cdptr;

  char XY[np][2] = { "P", "Q" };
  char X[np][2][2] = { {"A", "B"},
                       {"C", "D"} };

  for(int p=0; p<np; p++) {

    for(int braket=0; braket<2; braket++) {
      char XY_X_str[20];
      sprintf(XY_X_str,"%s%s",XY[p],X[p][braket]);
      rdptr vXY_X_ptr(new rdouble(XY_X_str));
      vXY_X[p][braket] = vXY_X_ptr;

      char W_XY_str[20];
      sprintf(W_XY_str,"W%s",XY[p]);
      rdptr vW_XY_ptr(new rdouble(W_XY_str));
      vW_XY[p] = vW_XY_ptr;

      char xyz_str[] = "xyz";
      for(int xyz=0; xyz<3; xyz++) {

        char XY_X_i_str[20];
        sprintf(XY_X_i_str,"%s_%c",XY_X_str,xyz_str[xyz]);
        rdptr XY_X_i_ptr(new rdouble(XY_X_i_str));
        XY_X[p][braket][xyz] = XY_X_i_ptr;

        char W_XY_i_str[20];
        sprintf(W_XY_i_str,"%s_%c",W_XY_str,xyz_str[xyz]);
        rdptr W_XY_i_ptr(new rdouble(W_XY_i_str));
        W_XY[p][xyz] = W_XY_i_ptr;

      }
    }

      char X_Y_str[20];
      sprintf(X_Y_str,"%s%s",X[p][0],X[p][1]);
      rdptr vX_Y_ptr(new rdouble(X_Y_str));
      vX_Y[p] = vX_Y_ptr;
      char xyz_str[] = "xyz";
      for(int xyz=0; xyz<3; xyz++) {
        char X_Y_i_str[20];
        sprintf(X_Y_i_str,"%s_%c",X_Y_str,xyz_str[xyz]);
        rdptr X_Y_i_ptr(new rdouble(X_Y_i_str));
        X_Y[p][xyz] = X_Y_i_ptr;
      }

  }
  
  //
  // TwoPRep ITR prefactors
  //

  for(int p=0; p<np; p++) {
    char tmp_str[40];
    sprintf(tmp_str,"TwoPRepITR_pfac0_%d",p);
    rdptr vpfac0_ptr(new rdouble(tmp_str));
    TwoPRepITR_vpfac0[p] = vpfac0_ptr;
    
    sprintf(tmp_str,"TwoPRepITR_pfac1_%d",p);
    rdptr pfac1_ptr(new rdouble(tmp_str));
    TwoPRepITR_pfac1[p] = pfac1_ptr;
    
    char xyz_str[] = "xyz";
    for(int xyz=0; xyz<3; xyz++) {
      char tmp_str[40];
      sprintf(tmp_str,"TwoPRepITR_pfac0_%d_%c",p,xyz_str[xyz]);
      rdptr pfac0_ptr(new rdouble(tmp_str));
      TwoPRepITR_pfac0[p][xyz] = pfac0_ptr;
    }
  }

  //
  // R12_k_G12 VRR prefactors
  //

  rdptr pfac2_ptr(new rdouble("R12kG12_pfac2"));
  R12kG12VRR_pfac2 = pfac2_ptr;

  for(int p=0; p<np; p++) {
    char tmp_str[20];
    sprintf(tmp_str,"R12kG12_pfac0_%d",p);
    rdptr vpfac0_ptr(new rdouble(tmp_str));
    R12kG12VRR_vpfac0[p] = vpfac0_ptr;
    
    sprintf(tmp_str,"R12kG12_pfac1_%d",p);
    rdptr pfac1_ptr(new rdouble(tmp_str));
    R12kG12VRR_pfac1[p] = pfac1_ptr;
    
    sprintf(tmp_str,"R12kG12_pfac3_%d",p);
    rdptr pfac3_ptr(new rdouble(tmp_str));
    R12kG12VRR_pfac3[p] = pfac3_ptr;
    
    sprintf(tmp_str,"R12kG12_pfac4_%d",p);
    rdptr vpfac4_ptr(new rdouble(tmp_str));
    R12kG12VRR_vpfac4[p] = vpfac4_ptr;
    
    char xyz_str[] = "xyz";
    for(int xyz=0; xyz<3; xyz++) {
      char tmp_str[20];
      sprintf(tmp_str,"R12kG12_pfac0_%d_%c",p,xyz_str[xyz]);
      rdptr pfac0_ptr(new rdouble(tmp_str));
      R12kG12VRR_pfac0[p][xyz] = pfac0_ptr;
      
      sprintf(tmp_str,"R12kG12_pfac4_%d_%c",p,xyz_str[xyz]);
      rdptr pfac4_ptr(new rdouble(tmp_str));
      R12kG12VRR_pfac4[p][xyz] = pfac4_ptr;
    }
  }

  for(int p=0; p<np; p++) {
    for(int braket=0; braket<2; braket++) {
      char tmpstr[200];
      sprintf(tmpstr,"zeta_%c",X[p][braket][0]);
      rdptr zptr(new rdouble(tmpstr));
      zeta[p][braket] = zptr;
    }
  }
  
  for(int p=0; p<np; p++) {
    for(int braket=0; braket<2; braket++) {
      char tmpstr[200];
      sprintf(tmpstr,"zeta_%c_2",X[p][braket][0]);
      rdptr zptr(new rdouble(tmpstr));
      zeta2[p][braket] = zptr;
    }
  }
  
  char alpha12_str[np][20] = { "zeta", "eta" };
  char alpha12_char[np+1] = "ze";
  for(int p=0; p<np; p++) {
    rdptr alpha12_ptr(new rdouble(alpha12_str[p]));
    alpha12[p] = alpha12_ptr;

    char oo2z[20];
    sprintf(oo2z,"oo2%c",alpha12_char[p]);
    rdptr oo2alpha12(new rdouble(oo2z));
    one_o_2alpha12[p] = oo2alpha12;

    char roz[20];
    sprintf(roz,"ro%c",alpha12_char[p]);
    rdptr roz_ptr(new rdouble(roz));
    rho_o_alpha12[p] = roz_ptr;
  }

  for(int i=0; i<NMAX; i++) {
    N_i[i] = prefactor::Scalar((double)i);
  }
}

Prefactors::~Prefactors()
{
}

SafePtr<Prefactors::cdouble>
Prefactors::Cdouble(double a)
{
  return prefactor::Scalar(a);
#if 0
  SafePtr<cdouble> tmp(new cdouble(a));
  typedef CTimeSingletons<double>::ManagerType ManagerType;
  const ManagerType::value_type& result = CTimeSingletons<double>::Manager()->find(tmp);
  return result.second;
#endif
}

namespace libint2 {
Prefactors prefactors;
};
