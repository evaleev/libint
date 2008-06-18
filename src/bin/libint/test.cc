
#include <boost/bind.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <rr.h>
#include <dg.h>
#include <dg.templ.h>
#include <integral.h>
#include <iter.h>
#include <policy_spec.h>
#include <intset_to_ints.h>
#include <strategy.h>
#include <buildtest.h>
#include <master.h>

using namespace std;
using namespace libint2;

namespace {
  int try_main (int argc, char* argv[]);
  void test0();
  void test1();
  void test2();
  void test3();
  void test4();
  void test5();
  void test_cgshell_iter(const CGShell& sh);

  template <class Callback>
    void RunTest(Callback test, const std::string& descr, std::ostream& os = std::cout);
  template <class Integral>
    void RunBuildTest(const typename Integral::BasisFunctionType& f1,
		      const typename Integral::BasisFunctionType& f2,
		      const typename Integral::BasisFunctionType& f3,
		      const typename Integral::BasisFunctionType& f4,
		      unsigned int size_to_unroll);
  template <class Integral>
    void RunBuildTest(const typename Integral::BasisFunctionType& f1,
		      const typename Integral::BasisFunctionType& f2,
		      const typename Integral::BasisFunctionType& f3,
		      const typename Integral::BasisFunctionType& f4,
		      unsigned int m,
		      unsigned int size_to_unroll);
};

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

namespace {

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
  SafePtr<CompilationParameters> cparams;

  int try_main (int argc, char* argv[])
  {
    LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
    taskmgr.add("test");
    
    // initialize cparams
    SafePtr<CompilationParameters> tmpcparams(new CompilationParameters);
    cparams = tmpcparams;
    cparams->max_am("test",2);

    // set default dims
    ImplicitDimensions::set_default_dims(cparams);
 
#if 1
    RunTest(test0,"iterators");
#endif
#if 1
    RunTest(test1,"memory managers");
#endif
#if 1
    RunTest(test2,"integrals types");
#endif
#if 1
    RunTest(test3,"recurrence relations");
#endif
#if 1
    RunTest(test4,"ERI build");
#endif

    return 0;
  }

  template <class Callback>
    void RunTest(Callback test, const std::string& descr, std::ostream& os) {
    const char hrule[] = "-----------------------------------------------------------------------------------------";
      os << hrule << endl;
      os << " Starting test: " << descr << endl;
      os << hrule << endl << endl;
      test();
      os << hrule << endl;
      os << " Finished test: " << descr << endl;
      os << hrule << endl << endl;
    }

  void
  test0()
  {
    SafePtr<TwoPRep_11_11_sq> pppp_quartet = TwoPRep_11_11_sq::Instance(sh_p,sh_p,sh_p,sh_p,0u);
    SafePtr<DGVertex> pppp_ptr = dynamic_pointer_cast<DGVertex,TwoPRep_11_11_sq>(pppp_quartet);

    // test CGShell iterator
    test_cgshell_iter(sh_s);
    test_cgshell_iter(sh_p);
    test_cgshell_iter(sh_d);
    test_cgshell_iter(sh_f);
    test_cgshell_iter(sh_g);
    test_cgshell_iter(sh_h);
    test_cgshell_iter(sh_i);
    test_cgshell_iter(sh_k);
    test_cgshell_iter(sh_l);

    // test IntegralSet iterator
    SafePtr<TwoPRep_11_11_sq> obj(pppp_quartet);
    cout << obj->description() << endl;
    SubIteratorBase< TwoPRep_11_11_sq > siter2(obj);
    for(siter2.init(); siter2; ++siter2)
      cout << siter2.elem()->description() << endl;
  }


  void
  test1()
  {
  
    SafePtr<Strategy> strat(new Strategy);
    SafePtr<Tactic> tactic(new FirstChoiceTactic<DummyRandomizePolicy>);
    SafePtr<TwoPRep_11_11_sq> xsxs_quartet = TwoPRep_11_11_sq::Instance(sh_f,sh_f,sh_f,sh_f,0u);
    cout << "Building " << xsxs_quartet->description() << endl;
    SafePtr<DGVertex> xsxs_ptr = dynamic_pointer_cast<DGVertex,TwoPRep_11_11_sq>(xsxs_quartet);

    SafePtr<MemoryManagerFactory> mmfactory(new MemoryManagerFactory);

    for(int m = 0; m < MemoryManagerFactory::ntypes; m++) {
      SafePtr<DirectedGraph> dg_xxxx3(new DirectedGraph);
      dg_xxxx3->append_target(xsxs_ptr);
      dg_xxxx3->apply(strat,tactic);
      dg_xxxx3->optimize_rr_out();
      dg_xxxx3->traverse();
      SafePtr<CodeContext> context(new CppCodeContext(cparams));
      SafePtr<MemoryManager> memman = mmfactory->memman(m);
      std::basic_ofstream<char> devnull("/dev/null");
      dg_xxxx3->generate_code(context,memman,ImplicitDimensions::default_dims(),SafePtr<CodeSymbols>(new CodeSymbols),xsxs_quartet->label(),devnull,devnull);
      cout << "Using " << mmfactory->label(m) << ": max memory used = " << memman->max_memory_used() << endl;
      dg_xxxx3->reset();
    }

    // Test BestFitMemoryFactory with tight_fit > 0
    const unsigned int tf_max = 6;
    for(int tf = 1; tf <= tf_max; tf++) {
      for(int ex=1; ex>=0; ex--) {
	SafePtr<DirectedGraph> dg_xxxx3(new DirectedGraph);
	dg_xxxx3->append_target(xsxs_ptr);
	dg_xxxx3->apply(strat,tactic);
	dg_xxxx3->optimize_rr_out();
	dg_xxxx3->traverse();
	SafePtr<CodeContext> context(new CppCodeContext(cparams));
	SafePtr<MemoryManager> memman(new BestFitMemoryManager(ex,tf));
	std::basic_ofstream<char> devnull("/dev/null");
	dg_xxxx3->generate_code(context,memman,ImplicitDimensions::default_dims(),SafePtr<CodeSymbols>(new CodeSymbols),xsxs_quartet->label(),devnull,devnull);
	cout << "Using BestFitMemoryManager(" << (ex == 1 ? "true" : "false")
	     << "," << tf << "): max memory used = " << memman->max_memory_used() << endl;
	dg_xxxx3->reset();
      }
    }
  }

  template <class Integral>
    void RunBuildTest(const typename Integral::BasisFunctionType& f1,
		      const typename Integral::BasisFunctionType& f2,
		      const typename Integral::BasisFunctionType& f3,
		      const typename Integral::BasisFunctionType& f4,
		      unsigned int size_to_unroll)
    {
      std::string descr("build ");
      SafePtr<Integral> i = Integral::Instance(f1,f2,f3,f4);
      descr += i->label();
      RunTest(boost::bind(__BuildTest<Integral,false>, i, cparams,
			  size_to_unroll, boost::ref(cout), SafePtr<Tactic>(new FirstChoiceTactic<DummyRandomizePolicy>),
			  SafePtr<MemoryManager>(new WorstFitMemoryManager), descr), descr);
    }

  template <class Integral>
    void RunBuildTest(const typename Integral::BasisFunctionType& f1,
		      const typename Integral::BasisFunctionType& f2,
		      const typename Integral::BasisFunctionType& f3,
		      const typename Integral::BasisFunctionType& f4,
		      unsigned int m, unsigned int size_to_unroll)
    {
      std::string descr("build ");
      SafePtr<Integral> i = Integral::Instance(f1,f2,f3,f4,m);
      descr += i->label();
      RunTest(boost::bind(__BuildTest<Integral,false>, i, cparams,
			  size_to_unroll, boost::ref(cout), SafePtr<Tactic>(new FirstChoiceTactic<DummyRandomizePolicy>),
			  SafePtr<MemoryManager>(new WorstFitMemoryManager), descr), descr);
    }


  void test_cgshell_iter(const CGShell& sh)
  {
    const unsigned int nbf = sh.num_bf();
    SafePtr<CGShell> sh_ptr(new CGShell(sh));
    sh_ptr->print(cout);
    SubIteratorBase<CGShell> siter1(sh_ptr);
    unsigned int bf = 0;
    for(siter1.init(); siter1; ++siter1, ++bf)
      siter1.elem().print(cout);
    if (bf != nbf)
      throw ProgrammingError("test::test_cgshell_iter -- number of basis functions from iterator and CGShell::num_bf do not match");
  }

  void
  test2()
  {
    {
      typedef TwoPRep_11_11_sq IType;
      SafePtr<IType> iset = IType::Instance(sh_p,sh_p,sh_p,sh_p,0u);
      std::cout << "Created integral set " << iset->label() << std::endl;
    }
    {
      typedef R12kG12_11_11_sq IType;
      SafePtr<IType> iset = IType::Instance(sh_p,sh_p,sh_p,sh_p,0u,IType::OperType(-1));
      std::cout << "Created integral set " << iset->label() << std::endl;
    }
  }

  void
  test3()
  {
    {
      typedef TwoPRep_11_11_sq IType;
      typedef VRR_a_11_TwoPRep_11_sh RRType;
      SafePtr<IType> iset = IType::Instance(sh_p,sh_p,sh_p,sh_p,0u);
      SafePtr<RRType> rr = RRType::Instance(iset,0);
      std::cout << "Created recurrence relation " << rr->label() << std::endl;
    }
    {
      typedef TwoPRep_11_11_sq IType;
      typedef VRR_c_11_TwoPRep_11_sh RRType;
      SafePtr<IType> iset = IType::Instance(sh_p,sh_p,sh_p,sh_p,0u);
      SafePtr<RRType> rr = RRType::Instance(iset,0);
      std::cout << "Created recurrence relation " << rr->label() << std::endl;
    }
  }

  void
  test4()
  {
    const unsigned int use_integrals = 1000000000;
    const unsigned int use_quartets = 1;
    RunBuildTest<TwoPRep_11_11_sq>(sh_p,sh_s,sh_p,sh_s,0,use_quartets);
    RunBuildTest<TwoPRep_11_11_sq>(sh_p,sh_p,sh_p,sh_p,0,use_quartets);
    RunBuildTest<TwoPRep_11_11_sq>(sh_p,sh_p,sh_p,sh_p,0,use_integrals);
  }
};
  
