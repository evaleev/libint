
#include <iostream>
#include <libint2.h>
#include <eri.h>
#include <prep_libint2.h>
#include <__ss_1_over_r_12_ds___up_0.h>

using namespace std;

int main(int argc, char** argv)
{

  typedef unsigned int uint;
  
  uint l1 = 0;
  uint m1 = 0;
  uint n1 = 0;
  uint l2 = 0;
  uint m2 = 0;
  uint n2 = 0;
  uint l3 = 2;
  uint m3 = 0;
  uint n3 = 0;
  uint l4 = 0;
  uint m4 = 0;
  uint n4 = 0;

  double alpha[4] = {0.5, 1.0, 1.5, 2.0};
  
  double A[3] = {1.0, 2.0, 3.0};
  double B[3] = {1.5, 2.5, 3.5};
  double C[3] = {4.0, 2.0, 0.0};
  double D[3] = {3.0, 3.0, 1.0};
  
  double ref_eri = eri(l1,m1,n1,alpha[0],A,
  l2,m2,n2,alpha[1],B,
  l3,m3,n3,alpha[2],C,
  l4,m4,n4,alpha[3],D,0);
  
  Libint_t* libint = new Libint_t;
  libint->stack = new double[LIBINT_MAX_STACK];
  prep_libint2(libint,0,alpha[0],A,
  0,alpha[1],B,
  0,alpha[2],C,
  0,alpha[3],D,0);

  cout << "Ref. eri = " << ref_eri << endl;
  
  //double new_eri = libint->__ss_1_over_r_12_ss___up_0;
  //compute__ss_1_over_r_12_ps___up_0(libint);
  //compute__ps_1_over_r_12_ps___up_0(libint);
  //compute__ds_1_over_r_12_ds___up_0(libint);
  compute__ss_1_over_r_12_ds___up_0(libint);
  //compute__ps_1_over_r_12_ds___up_0(libint);
  //compute__ds_1_over_r_12_ps___up_0(libint);
  //compute__ds_1_over_r_12_ss___up_0(libint);
  //compute__fs_1_over_r_12_ss___up_0(libint);
  //compute__ss_1_over_r_12_fs___up_0(libint);
  for(int i=0;i<6; i++)
    cout << "New  eri[" << i << "] = " << libint->targets[0][i] << endl;
}


