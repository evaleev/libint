
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <math.h>

#include <rr.h>

#include <libint2.h>
#include <prep_libint2v.h>
#include <__ds_1_over_r_12_ds___up_0.h>

using namespace std;
using namespace libint2;

namespace {
  std::string usage();
};

int main(int argc, char** argv)
{
  typedef unsigned int uint;

  if (argc != 3) {
    cerr << usage() << endl;
    exit(1);
  }

  uint niter = atoi(argv[1]);
  uint veclength = atoi(argv[2]);

  uint am[4] = {2, 0, 2, 0};
  double alpha[4] = {0.5, 1.0, 1.5, 2.0};
  double A[3] = {1.0, 2.0, 3.0};
  double B[3] = {1.5, 2.5, 3.5};
  double C[3] = {4.0, 2.0, 0.0};
  double D[3] = {3.0, 3.0, 1.0};

  SafePtr<CGShell> sh0(new CGShell(&(am[0])));
  SafePtr<CGShell> sh1(new CGShell(&(am[1])));
  SafePtr<CGShell> sh2(new CGShell(&(am[2])));
  SafePtr<CGShell> sh3(new CGShell(&(am[3])));
  
  Libint_t* libint = new Libint_t;
  libint->stack = new double[10000];
  prep_libint2v(libint,veclength,am[0],alpha[0],A,
  am[1],alpha[1],B,
  am[2],alpha[2],C,
  am[3],alpha[3],D,0);

  cout << "Computing (" << sh0->label() << sh1->label()
       << "|" << sh2->label() << sh3->label() << ") " << niter << " times" << endl;

  for(int iter=0; iter<niter; iter++)
    compute__ds_1_over_r_12_ds___up_0(libint);

  exit(0);
}

namespace {
  std::string
  usage()
  {
    ostringstream oss;
    oss << "Usage: time_eri_vector <niter> <vecl>" << endl;
    oss << "       niter -- number of iterations" << endl;
    oss << "       vecl  -- vector length" << endl;
    return oss.str();
  }
};

