
#include <iostream>
#include <vector>
#include <rr.h>
#include <dg.h>

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
  //cout << "Permutational symmetry of 2e repulsion = " << TwoERep.psymm(0,1) << endl;
  cerr << "ok" << endl;

  cerr << "Testing GaussShell ... ";
  unsigned int a[] = {1, 0, 0};
  GaussShell A(a);
  cerr << "ok" << endl;

  // Need to be able to compile these constructors

  // Make a (ds|ps) = <dp|ss> class
  unsigned int am[][1] = { {0}, {1}, {2}, {3} };
  CGShell sh_s(am[0]);
  CGShell sh_p(am[1]);
  CGShell sh_d(am[2]);
  CGShell sh_f(am[3]);
  vector<CGShell> sh_12[TwoERep::np];
  sh_12[0].push_back(sh_d);
  sh_12[1].push_back(sh_p);
  vector<CGShell> sh_34[TwoERep::np];
  sh_34[0].push_back(sh_s);
  sh_34[1].push_back(sh_s);

  TwoERep_2b2k<CGShell> sq_a(sh_12,sh_34, 0);
  sq_a.print();
  VRR_a_ERI_2b2k_shell vrr1(&sq_a);

  TwoERep_2b2k<CGShell> sq_b(sh_34,sh_34, 0);
  sq_b.print();
  VRR_c_ERI_2b2k_shell vrr2(&sq_b);

  // Create a DAG for a VRR case
  DirectedGraph dg_vrr_dsps;
  dg_vrr_dsps.append_target< TwoERep_2b2k<CGShell>, VRR_a_ERI_2b2k_shell>(&sq_a);

}

