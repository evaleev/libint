
#include <iostream>
#include <vector>
#include <rr.h>
#include <dg.h>
#include <typelist.h>
#include <integral.h>

using namespace std;
using namespace libint2;

static int try_main (int argc, char* argv[]);
static void test_typelists();
static void test_old_integrals();
static void test_integrals();

int main (int argc, char* argv[])
{

  try {
    try_main(argc,argv);
  }
  catch(std::exception& a) {
    cout << endl
         << "  WARNING! Caught a standard exception:" << endl
         << "    " << a.what() << endl << endl;
  }
}

int try_main (int argc, char* argv[])
{

  unsigned int am[][1] = { {0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}};
  CGShell sh_s(am[0]);
  CGShell sh_p(am[1]);
  CGShell sh_d(am[2]);
  CGShell sh_f(am[3]);
  CGShell sh_g(am[4]);
  CGShell sh_h(am[5]);
  CGShell sh_i(am[6]);
  CGShell sh_k(am[7]);
  CGShell sh_l(am[8]);
  CGShell sh_m(am[9]);
  CGShell sh_n(am[10]);

  vector<CGShell> sh_12[TwoERep::np];
  sh_12[0].push_back(sh_p);
  sh_12[1].push_back(sh_p);
  vector<CGShell> sh_34[TwoERep::np];
  sh_34[0].push_back(sh_p);
  sh_34[1].push_back(sh_p);

  TwoERep_2b2k<CGShell>* sq1 = TwoERep_2b2k<CGShell>::Instance(sh_12,sh_34, 0);
  sq1->print();
  VRR_c_ERI_2b2k_shell vrr1(sq1);

  // Create a DAG for a VRR case
  DirectedGraph dg_vrr_xxxx;

  dg_vrr_xxxx.append_target(sq1);
  dg_vrr_xxxx.apply_to_all< HRR_ab_ERI_2b2k_shell >();
  dg_vrr_xxxx.apply_to_all< HRR_cd_ERI_2b2k_shell >();
  dg_vrr_xxxx.apply_to_all< VRR_a_ERI_2b2k_shell >();
  dg_vrr_xxxx.apply_to_all< VRR_c_ERI_2b2k_shell >();

  dg_vrr_xxxx.traverse();
  dg_vrr_xxxx.debug_print_traversal(cout);

  test_typelists();
  
  typedef TwoPRep_11_11<CGShell> TwoPRep_sh_11_11;
  TwoPRep_sh_11_11* pppp_quartet = TwoPRep_sh_11_11::Instance(sh_p,sh_p,sh_p,sh_p,0);
  pppp_quartet->print(cout);

  // Create a DAG for a VRR case
  DirectedGraph dg_xxxx;

  dg_xxxx.append_target(pppp_quartet);
  dg_xxxx.apply_to_all< HRR_ba_11_TwoPRep_11_sh >();
  dg_xxxx.apply_to_all< HRR_dc_11_TwoPRep_11_sh >();
  dg_xxxx.apply_to_all< VRR_a_11_TwoPRep_11_sh >();
  dg_xxxx.apply_to_all< VRR_c_11_TwoPRep_11_sh >();
  dg_xxxx.apply_to_all< VRR_b_11_TwoPRep_11_sh >();
  dg_xxxx.apply_to_all< VRR_d_11_TwoPRep_11_sh >();

  dg_xxxx.traverse();
  dg_xxxx.debug_print_traversal(cout);


}


void
test_typelists()
{
  typedef PTYPELIST_2( Tag2, CGShell, CGShell) two_shells;
  two_shells xx;
  typedef PTYPELIST_3( Tag2, CGShell, CGShell, CGShell) three_shells;
  three_shells xxx;
  typedef GenIntegral<TwoERep,two_shells,two_shells> NewTwoERep_2b2k;
  NewTwoERep_2b2k xxxx(xx,xx);
  
  ValueAt<two_shells, 1>::ResultType xx_1 = ValueAt<two_shells,1>::typelist_member(xx);

  cout << "Length of my_typelist = " << Length<two_shells>::value << endl;
  CGShell xx_2(xx_1);
  xx_2.inc();
}

