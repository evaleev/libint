
#include <smart_ptr.h>
#include <prefactors.h>

using namespace libint2;

Prefactors::Prefactors() :
  rho(new rdouble("rho")),
  one_o_2alphasum(new rdouble("oo2ze"))
{
  typedef SafePtr<rdouble> rdptr;
  typedef SafePtr<rdouble> cdptr;

  char XY[np][3] = { "P", "Q" };
  char X[np][2][3] = { {"A", "B"},
                       {"C", "D"} };

  for(int p=0; p<np; p++)
    for(int braket=0; braket<1; braket++) {
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
    sprintf(roz,"ro2%c",alpha12_char[p]);
    rdptr roz_ptr(new rdouble(roz));
    rho_o_alpha12[p] = roz_ptr;
  }

  for(int i=0; i<NMAX; i++) {
    double i_fp = (double) i;
    char c_str_repr[20];
    sprintf(c_str_repr,"%.2lf", i_fp);
    SafePtr<cdouble> N_ptr(new cdouble(c_str_repr,i_fp));
    N_i[i] = N_ptr;
  }
}

Prefactors::~Prefactors()
{
}

namespace libint2 {
Prefactors prefactors;
};
