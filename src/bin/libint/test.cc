
#include <iostream>
#include <rr.h>

using namespace std;
using namespace libint2;

int main (int argc, char* argv[])
{

  cerr << "Testing CartesianGaussian ... ";
  CartesianGaussian p1(0, 0, 0);
  CartesianGaussian p2 = p1;
  cerr << "ok" << endl;

  cerr << "Testing Operator ... ";
  Operator OneENucAttr("One-Electron Nuclear Attraction","1ENA",1,vector<char>(0,0));
  cout << "Permutational symmetry of 2e repulsion = " << TwoERep.psymm(0,1) << endl;
  cerr << "ok" << endl;

  cerr << "Testing GaussShell ... ";
  unsigned int a[] = {1, 0, 0};
  GaussShell A(a);
  cerr << "ok" << endl;
  

}

