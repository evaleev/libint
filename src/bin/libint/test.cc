
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
#include <master_ints_list.h>
#include <master_rrs_list.h>

using namespace std;
using namespace libint2;

long living_count = 0;

namespace {
  int try_main (int argc, char* argv[]);
  void test0();
  void test1();
  void test2();
  void test3();
  void test4();
  void test5();
  void test6();
  void test7();
  void test8();
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
  template <class Integral>
    void RunBuildTest(const typename Integral::BasisFunctionType& f1,
              const typename Integral::BasisFunctionType& f2,
              const typename Integral::BasisFunctionType& f3,
              const typename Integral::BasisFunctionType& f4,
              unsigned int m,
              const typename Integral::OperType::Descriptor& descr,
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
  return 0;
}

namespace {

  CGShell sh_s(0);
  CGShell sh_p(1);
  CGShell sh_d(2);
  CGShell sh_f(3);
  CGShell sh_g(4);
  CGShell sh_h(5);
  CGShell sh_i(6);
  CGShell sh_k(7);
  CGShell sh_l(8);
  CGShell sh_m(9);
  CGShell sh_n(10);
  CGShell sh_o(11);
  CGShell sh_q(12);
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

#if 0
    RunTest(test0,"iterators");
#endif
#if 0
    RunTest(test1,"memory managers");
#endif
#if 0
    RunTest(test2,"integrals types");
#endif
#if 0
    RunTest(test3,"recurrence relations");
#endif
#if 0
    RunTest(test4,"primitive ERI build");
#endif
#if 0
    RunTest(test5,"contracted ERI build");
#endif
#if 0
    RunTest(test6,"contracted derivative ERI build");
#endif
#if 0
    RunTest(test7,"shell-set RR generation");
#endif
#if 1
    RunTest(test8,"congracted G12 integral build");
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

    for(unsigned int m = 0; m < MemoryManagerFactory::ntypes; m++) {
      SafePtr<DirectedGraph> dg_xxxx3(new DirectedGraph);
      SafePtr<CodeContext> context(new CppCodeContext(cparams));
      SafePtr<MemoryManager> memman = mmfactory->memman(m);
      dg_xxxx3->append_target(xsxs_ptr);
      dg_xxxx3->apply(strat,tactic);
      dg_xxxx3->optimize_rr_out(context);
      dg_xxxx3->traverse();
      std::basic_ofstream<char> devnull("/dev/null");
      dg_xxxx3->generate_code(context,memman,ImplicitDimensions::default_dims(),SafePtr<CodeSymbols>(new CodeSymbols),xsxs_quartet->label(),devnull,devnull);
      cout << "Using " << mmfactory->label(m) << ": max memory used = " << memman->max_memory_used() << endl;
      dg_xxxx3->reset();
    }

    // Test BestFitMemoryFactory with tight_fit > 0
    const unsigned int tf_max = 6;
    for(unsigned int tf = 1; tf <= tf_max; tf++) {
      for(int ex=1; ex>=0; ex--) {
	SafePtr<DirectedGraph> dg_xxxx3(new DirectedGraph);
    SafePtr<CodeContext> context(new CppCodeContext(cparams));
    SafePtr<MemoryManager> memman(new BestFitMemoryManager(ex,tf));
	dg_xxxx3->append_target(xsxs_ptr);
	dg_xxxx3->apply(strat,tactic);
	dg_xxxx3->optimize_rr_out(context);
	dg_xxxx3->traverse();
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
      std::vector<SafePtr<Integral> > targets(1,i);
      RunTest(boost::bind(__BuildTest<Integral,false>, targets, cparams,
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
      std::vector<SafePtr<Integral> > targets(1,i);
      RunTest(boost::bind(__BuildTest<Integral,false>, targets, cparams,
			  size_to_unroll, boost::ref(cout), SafePtr<Tactic>(new FirstChoiceTactic<DummyRandomizePolicy>),
			  SafePtr<MemoryManager>(new WorstFitMemoryManager), descr), descr);
    }
  template <class Integral>
    void RunBuildTest(const typename Integral::BasisFunctionType& f1,
              const typename Integral::BasisFunctionType& f2,
              const typename Integral::BasisFunctionType& f3,
              const typename Integral::BasisFunctionType& f4,
              unsigned int m,
              const typename Integral::OperType::Descriptor& descr, unsigned int size_to_unroll)
    {
      std::string descr_label("build ");
      typedef typename Integral::OperType::Descriptor Descriptor;
      GenOper<Descriptor> oper(descr);
      SafePtr<Integral> i = Integral::Instance(f1,f2,f3,f4,m,oper);
      descr_label += i->label();
      std::vector<SafePtr<Integral> > targets(1,i);
      RunTest(boost::bind(__BuildTest<Integral,false>, targets, cparams,
              size_to_unroll, boost::ref(cout), SafePtr<Tactic>(new FirstChoiceTactic<DummyRandomizePolicy>),
              SafePtr<MemoryManager>(new WorstFitMemoryManager), descr_label), descr_label);
    }


  void test_cgshell_iter(const CGShell& sh)
  {
    const unsigned int nbf = sh.num_bf();
    SafePtr<CGShell> sh_ptr(new CGShell(sh));
    sh_ptr->print(cout);
    SubIteratorBase<CGShell> siter1(*sh_ptr);
    unsigned int bf = 0;
    for(siter1.init(); siter1; ++siter1, ++bf)
      siter1.elem().print(cout);
    if (bf != nbf)
      throw ProgrammingError("test::test_cgshell_iter -- number of basis functions from iterator and CGShell::num_bf do not match");
  }

  void
  test2()
  {
    CGShell::set_contracted_default_value(true);
    const bool contracted_targets_old_value = cparams->contracted_targets();
    cparams->contracted_targets(true);
    CGShell csh_s(0u);
    CGShell csh_p(1u);
    CGShell csh_d(2u);
    CGShell csh_d_dx(2u); csh_d_dx.deriv().inc(0,1);
    CGShell csh_f_dx(3u); csh_f_dx.deriv().inc(0,1);
    CGShell csh_q(12u);

    {
      typedef TwoPRep_11_11_sq IType;
      SafePtr<IType> iset = IType::Instance(sh_p,sh_p,sh_p,sh_p,0u);
      std::cout << "Created integral set " << iset->label() << " key = " << iset->key() << std::endl;
    }
    {
      typedef R12kG12_11_11_sq IType;
      SafePtr<IType> iset = IType::Instance(sh_p,sh_p,sh_p,sh_p,0u,IType::OperType(-1));
      std::cout << "Created integral set " << iset->label() << " key = " << iset->key() << std::endl;
    }
    {
      typedef TwoPRep_11_11_sq IType;
      SafePtr<IType> iset = IType::Instance(csh_s,csh_q,csh_s,csh_s,0u);
      std::cout << "Created integral set " << iset->label() << " key = " << iset->key() << std::endl;
    }
    {
      typedef TwoPRep_11_11_sq IType;
      SafePtr<IType> iset = IType::Instance(csh_s,csh_d_dx,csh_s,csh_s,0u);
      std::cout << "Created integral set " << iset->label() << " key = " << iset->key() << std::endl;
    }
    {
      typedef TwoPRep_11_11_sq IType;
      SafePtr<IType> iset = IType::Instance(csh_s,csh_f_dx,csh_s,csh_s,0u);
      std::cout << "Created integral set " << iset->label() << " key = " << iset->key() << std::endl;
    }
    {
      typedef TwoPRep_11_11_sq IType;
      SafePtr<IType> iset = IType::Instance(csh_q,csh_q,csh_q,csh_q,0u);
      std::cout << "Created integral set " << iset->label() << " key = " << iset->key() << std::endl;
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
    {
      typedef DivG12prime_xTx_11_11_sq IType;
      typedef CR_11_DivG12prime_xTx_11_sh RRType;
      SafePtr<IType> iset = IType::Instance(sh_p,sh_p,sh_p,sh_p,0u,DivG12prime_xTx_Descr(0));
      SafePtr<RRType> rr = RRType::Instance(iset,0);
      std::cout << "Created recurrence relation " << rr->label() << std::endl;
    }
  }

  // primitive ERI build
  void
  test4()
  {
    const unsigned int use_integrals = 1000000000;
    const unsigned int use_quartets = 1;
    RunBuildTest<TwoPRep_11_11_sq>(sh_p,sh_s,sh_p,sh_s,0,use_quartets);
    RunBuildTest<TwoPRep_11_11_sq>(sh_p,sh_p,sh_p,sh_p,0,use_quartets);
    RunBuildTest<TwoPRep_11_11_sq>(sh_p,sh_p,sh_p,sh_p,0,use_integrals);
  }

  // contracted ERI build
  void
  test5()
  {
    const unsigned int use_integrals = 1000000000;
    const unsigned int use_quartets = 1;

    CGShell::set_contracted_default_value(true);
    const bool contracted_targets_old_value = cparams->contracted_targets();
    cparams->contracted_targets(true);
    CGShell csh_s(0u);
    CGShell csh_p(1u);

    //RunBuildTest<TwoPRep_11_11_sq>(csh_p,csh_s,csh_p,csh_s,0,use_quartets);
    RunBuildTest<TwoPRep_11_11_sq>(csh_s,csh_p,csh_s,csh_s,0,use_quartets);

    cparams->contracted_targets(contracted_targets_old_value);
  }

  // contracted derivative ERI build
  void
  test6()
  {
    const unsigned int use_integrals = 1000000000;
    const unsigned int use_quartets = 1;

    CGShell::set_contracted_default_value(true);
    const bool contracted_targets_old_value = cparams->contracted_targets();
    cparams->contracted_targets(true);
    CGShell csh_s(0u);
    CGShell csh_p(1u);
    CGShell csh_d(2u);
    CGShell csh_s_dx(0u); csh_s_dx.deriv().inc(0,1);
    CGShell csh_p_dx(1u); csh_p_dx.deriv().inc(0,1);
    CGShell csh_p_dy(1u); csh_p_dy.deriv().inc(1,1);
    CGShell csh_p_dxyz(1u); csh_p_dxyz.deriv().inc(0,1); csh_p_dxyz.deriv().inc(1,1); csh_p_dxyz.deriv().inc(2,1);
    CGShell csh_d_dx(2u); csh_d_dx.deriv().inc(0,1);
    CGShell csh_s_d2x(0u); csh_s_d2x.deriv().inc(0,2);
    CGShell csh_p_d2x(1u); csh_p_d2x.deriv().inc(0,2);

    RunBuildTest<TwoPRep_11_11_sq>(csh_d_dx,csh_s,csh_s,csh_s,0,use_integrals);
    //RunBuildTest<TwoPRep_11_11_sq>(csh_d_dx,csh_s,csh_s,csh_s,0,use_quartets);

    cparams->contracted_targets(contracted_targets_old_value);
  }

  // testing RR generators
  void
  test7()
  {
    CGShell::set_contracted_default_value(true);
    const bool contracted_targets_old_value = cparams->contracted_targets();
    cparams->contracted_targets(true);

    CGShell csh_s(0u);
    CGShell csh_p(1u);
    CGShell csh_d(2u);
    CGShell csh_s_dx(0u); csh_s_dx.deriv().inc(0,1);
    CGShell csh_p_dx(1u); csh_p_dx.deriv().inc(0,1);
    CGShell csh_p_dy(1u); csh_p_dy.deriv().inc(1,1);
    CGShell csh_p_dxyz(1u); csh_p_dxyz.deriv().inc(0,1); csh_p_dxyz.deriv().inc(1,1); csh_p_dxyz.deriv().inc(2,1);
    CGShell csh_d_dx(2u); csh_d_dx.deriv().inc(0,1);
    CGShell csh_s_d2x(0u); csh_s_d2x.deriv().inc(0,2);
    CGShell csh_p_d2x(1u); csh_p_d2x.deriv().inc(0,2);

    SafePtr<TwoPRep_11_11_sq> target = TwoPRep_11_11_sq::Instance(csh_p,csh_d_dx,csh_s,csh_s,mType(0u));

    SafePtr<RecurrenceRelation> rr = HRR_ba_11_TwoPRep_11_sh::Instance(target,0);
    assert(rr != 0 && rr->num_children() != 0);
    SafePtr<RRStack> rrstack = RRStack::Instance();
    rrstack->find(rr);

    std::deque<std::string> decl_filenames, def_filenames;
    generate_rr_code(std::cout,
                     cparams,
                     decl_filenames, def_filenames);

    cparams->contracted_targets(contracted_targets_old_value);
  }

  // contracted g12 integral build
  void
  test8()
  {
    const unsigned int use_integrals = 1000000000;
    const unsigned int use_quartets = 1;

    CGShell::set_contracted_default_value(true);
    const bool contracted_targets_old_value = cparams->contracted_targets();
    cparams->contracted_targets(true);
    CGShell csh_s(0u);
    CGShell csh_p(1u);

    const Ti_G12_Descr t0g12_descr(0);
    RunBuildTest<TiG12_11_11_sq>(csh_s,csh_s,csh_s,csh_s,0,t0g12_descr,use_quartets);

    cparams->contracted_targets(contracted_targets_old_value);
  }

};

