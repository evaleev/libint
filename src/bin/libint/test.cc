/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint compiler.
 *
 *  Libint compiler is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint compiler is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint compiler.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <buildtest.h>
#include <dg.h>
#include <integral.h>
#include <intset_to_ints.h>
#include <iter.h>
#include <master_ints_list.h>
#include <master_rrs_list.h>
#include <policy_spec.h>
#include <rr.h>
#include <strategy.h>

#include <boost/bind.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace libint2;

long living_count = 0;
const unsigned int use_integrals = std::numeric_limits<unsigned int>::max();
const unsigned int use_quartets = 0;

namespace {
int try_main(int argc, char* argv[]);
void test0();
void test1();
void test2();
void test3();
void test4();
void test5();
void test6();
void test7();
void test8();
void test9();
void test_cgshell_iter(const CGShell& sh);

template <class Callback>
void RunTest(Callback test, const std::string& descr,
             std::ostream& os = std::cout);
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
                  unsigned int m, unsigned int size_to_unroll);
template <class Integral>
void RunBuildTest(const typename Integral::BasisFunctionType& f1,
                  const typename Integral::BasisFunctionType& f2,
                  const typename Integral::BasisFunctionType& f3,
                  const typename Integral::BasisFunctionType& f4,
                  unsigned int m,
                  const typename Integral::OperType::Descriptor& descr,
                  unsigned int size_to_unroll);
};  // namespace

int main(int argc, char* argv[]) {
  try {
    try_main(argc, argv);
  } catch (std::exception& a) {
    cout << endl
         << "  WARNING! Caught a standard exception:" << endl
         << "    " << a.what() << endl
         << endl;
    return 1;
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
std::shared_ptr<CompilationParameters> cparams;

int try_main(int argc, char* argv[]) {
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.add("test");

  // initialize cparams
  std::shared_ptr<CompilationParameters> tmpcparams(new CompilationParameters);
  cparams = tmpcparams;
  cparams->max_am("test", 2);

  // set default dims
  ImplicitDimensions::set_default_dims(cparams);

#if 1
  RunTest(test0, "iterators");
#endif
#if 1
  RunTest(test1, "memory managers");
#endif
#if 1
  RunTest(test2, "integrals types");
#endif
#if 1
  RunTest(test3, "recurrence relations");
#endif
#if 1
  RunTest(test4, "primitive ERI build");
#endif
#if 1
  RunTest(test5, "contracted ERI build");
#endif
#if 1
  RunTest(test6, "contracted derivative ERI build");
#endif
#if 1
  RunTest(test7, "shell-set RR generation");
#endif
#if 1
  RunTest(test8, "contracted G12 integral build");
#endif
#if 1
  RunTest(test9, "DerivGaussV2 relation");
#endif

  return 0;
}

template <class Callback>
void RunTest(Callback test, const std::string& descr, std::ostream& os) {
  const char hrule[] =
      "------------------------------------------------------------------------"
      "-----------------";
  os << hrule << endl;
  os << " Starting test: " << descr << endl;
  os << hrule << endl << endl;
  test();
  os << hrule << endl;
  os << " Finished test: " << descr << endl;
  os << hrule << endl << endl;
}

void test0() {
  // test CGShell labels
  if (CGShell(0).label() != "s")
    throw ProgrammingError("CGShell::label() failed for l=0");
  if (CGShell(20).label() != "z")
    throw ProgrammingError("CGShell::label() failed for l=20");
  if (CGShell(21).label() != "ps")
    throw ProgrammingError("CGShell::label() failed for l=21");
  if (CGShell(22).label() != "pp")
    throw ProgrammingError("CGShell::label() failed for l=22");
  if (CGShell(42).label() != "ds")
    throw ProgrammingError("CGShell::label() failed for l=42");

  std::shared_ptr<TwoPRep_11_11_sq> pppp_quartet =
      TwoPRep_11_11_sq::Instance(sh_p, sh_p, sh_p, sh_p, 0u);
  std::shared_ptr<DGVertex> pppp_ptr =
      std::dynamic_pointer_cast<DGVertex, TwoPRep_11_11_sq>(pppp_quartet);

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
  std::shared_ptr<TwoPRep_11_11_sq> obj(pppp_quartet);
  cout << obj->description() << endl;
  SubIteratorBase<TwoPRep_11_11_sq> siter2(obj);
  for (siter2.init(); siter2; ++siter2)
    cout << siter2.elem()->description() << endl;
}

void test1() {
  std::shared_ptr<Strategy> strat(new Strategy);
  std::shared_ptr<Tactic> tactic(new FirstChoiceTactic<DummyRandomizePolicy>);
  std::shared_ptr<TwoPRep_11_11_sq> xsxs_quartet =
      TwoPRep_11_11_sq::Instance(sh_f, sh_f, sh_f, sh_f, 0u);
  cout << "Building " << xsxs_quartet->description() << endl;
  std::shared_ptr<DGVertex> xsxs_ptr =
      std::dynamic_pointer_cast<DGVertex, TwoPRep_11_11_sq>(xsxs_quartet);

  std::shared_ptr<MemoryManagerFactory> mmfactory(new MemoryManagerFactory);

  for (unsigned int m = 0; m < MemoryManagerFactory::ntypes; m++) {
    std::shared_ptr<DirectedGraph> dg_xxxx3(new DirectedGraph);
    std::shared_ptr<CodeContext> context(new CppCodeContext(cparams));
    std::shared_ptr<MemoryManager> memman = mmfactory->memman(m);
    dg_xxxx3->append_target(xsxs_ptr);
    dg_xxxx3->apply(strat, tactic);
    dg_xxxx3->optimize_rr_out(context);
    dg_xxxx3->traverse();
    std::basic_ofstream<char> devnull("/dev/null");
    dg_xxxx3->generate_code(context, memman, ImplicitDimensions::default_dims(),
                            std::shared_ptr<CodeSymbols>(new CodeSymbols),
                            xsxs_quartet->label(), devnull, devnull);
    cout << "Using " << mmfactory->label(m)
         << ": max memory used = " << memman->max_memory_used() << endl;
    dg_xxxx3->reset();
  }

  // Test BestFitMemoryFactory with tight_fit > 0
  const unsigned int tf_max = 6;
  for (unsigned int tf = 1; tf <= tf_max; tf++) {
    for (int ex = 1; ex >= 0; ex--) {
      std::shared_ptr<DirectedGraph> dg_xxxx3(new DirectedGraph);
      std::shared_ptr<CodeContext> context(new CppCodeContext(cparams));
      std::shared_ptr<MemoryManager> memman(new BestFitMemoryManager(ex, tf));
      dg_xxxx3->append_target(xsxs_ptr);
      dg_xxxx3->apply(strat, tactic);
      dg_xxxx3->optimize_rr_out(context);
      dg_xxxx3->traverse();
      std::basic_ofstream<char> devnull("/dev/null");
      dg_xxxx3->generate_code(context, memman,
                              ImplicitDimensions::default_dims(),
                              std::shared_ptr<CodeSymbols>(new CodeSymbols),
                              xsxs_quartet->label(), devnull, devnull);
      cout << "Using BestFitMemoryManager(" << (ex == 1 ? "true" : "false")
           << "," << tf << "): max memory used = " << memman->max_memory_used()
           << endl;
      dg_xxxx3->reset();
    }
  }
}

template <class Integral>
void RunBuildTest(const typename Integral::BasisFunctionType& f1,
                  const typename Integral::BasisFunctionType& f2,
                  const typename Integral::BasisFunctionType& f3,
                  const typename Integral::BasisFunctionType& f4,
                  unsigned int size_to_unroll) {
  std::string descr("build ");
  std::shared_ptr<Integral> i = Integral::Instance(f1, f2, f3, f4);
  descr += i->label();
  std::vector<std::shared_ptr<Integral> > targets(1, i);
  RunTest(
      boost::bind(
          __BuildTest<Integral, false>, targets, cparams, size_to_unroll,
          boost::ref(cout),
          std::shared_ptr<Tactic>(new FirstChoiceTactic<DummyRandomizePolicy>),
          std::shared_ptr<MemoryManager>(new WorstFitMemoryManager), descr),
      descr);
}

template <class Integral>
void RunBuildTest(const typename Integral::BasisFunctionType& f1,
                  const typename Integral::BasisFunctionType& f2,
                  const typename Integral::BasisFunctionType& f3,
                  const typename Integral::BasisFunctionType& f4,
                  unsigned int m, unsigned int size_to_unroll) {
  std::string descr("build ");
  std::shared_ptr<Integral> i = Integral::Instance(f1, f2, f3, f4, m);
  descr += i->label();
  std::vector<std::shared_ptr<Integral> > targets(1, i);
  RunTest(
      boost::bind(
          __BuildTest<Integral, false>, targets, cparams, size_to_unroll,
          boost::ref(cout),
          std::shared_ptr<Tactic>(new FirstChoiceTactic<DummyRandomizePolicy>),
          std::shared_ptr<MemoryManager>(new WorstFitMemoryManager), descr),
      descr);
}
template <class Integral>
void RunBuildTest(const typename Integral::BasisFunctionType& f1,
                  const typename Integral::BasisFunctionType& f2,
                  const typename Integral::BasisFunctionType& f3,
                  const typename Integral::BasisFunctionType& f4,
                  unsigned int m,
                  const typename Integral::OperType::Descriptor& descr,
                  unsigned int size_to_unroll) {
  std::string descr_label("build ");
  typedef typename Integral::OperType::Descriptor Descriptor;
  GenOper<Descriptor> oper(descr);
  std::shared_ptr<Integral> i = Integral::Instance(f1, f2, f3, f4, m, oper);
  descr_label += i->label();
  std::vector<std::shared_ptr<Integral> > targets(1, i);
  RunTest(boost::bind(__BuildTest<Integral, false>, targets, cparams,
                      size_to_unroll, boost::ref(cout),
                      std::shared_ptr<Tactic>(
                          new FirstChoiceTactic<DummyRandomizePolicy>),
                      std::shared_ptr<MemoryManager>(new WorstFitMemoryManager),
                      descr_label),
          descr_label);
}

void test_cgshell_iter(const CGShell& sh) {
  const unsigned int nbf = sh.num_bf();
  std::shared_ptr<CGShell> sh_ptr(new CGShell(sh));
  sh_ptr->print(cout);
  SubIteratorBase<CGShell> siter1(*sh_ptr);
  unsigned int bf = 0;
  for (siter1.init(); siter1; ++siter1, ++bf) siter1.elem().print(cout);
  if (bf != nbf)
    throw ProgrammingError(
        "test::test_cgshell_iter -- number of basis functions from iterator "
        "and CGShell::num_bf do not match");
}

void test2() {
  CGShell::set_contracted_default_value(true);
  cparams->contracted_targets(true);
  CGShell csh_s(0u);
  CGShell csh_p(1u);
  CGShell csh_d(2u);
  CGShell csh_d_dx(2u);
  csh_d_dx.deriv().inc(0, 1);
  CGShell csh_f_dx(3u);
  csh_f_dx.deriv().inc(0, 1);
  CGShell csh_q(12u);

  {
    typedef TwoPRep_11_11_sq IType;
    std::shared_ptr<IType> iset = IType::Instance(sh_p, sh_p, sh_p, sh_p, 0u);
    std::cout << "Created integral set " << iset->label()
              << " key = " << iset->key() << std::endl;
  }
  {
    typedef R12kG12_11_11_sq IType;
    std::shared_ptr<IType> iset =
        IType::Instance(sh_p, sh_p, sh_p, sh_p, 0u, IType::OperType(-1));
    std::cout << "Created integral set " << iset->label()
              << " key = " << iset->key() << std::endl;
  }
  {
    typedef TwoPRep_11_11_sq IType;
    std::shared_ptr<IType> iset =
        IType::Instance(csh_s, csh_q, csh_s, csh_s, 0u);
    std::cout << "Created integral set " << iset->label()
              << " key = " << iset->key() << std::endl;
  }
  {
    typedef TwoPRep_11_11_sq IType;
    std::shared_ptr<IType> iset =
        IType::Instance(csh_s, csh_d_dx, csh_s, csh_s, 0u);
    std::cout << "Created integral set " << iset->label()
              << " key = " << iset->key() << std::endl;
  }
  {
    typedef TwoPRep_11_11_sq IType;
    std::shared_ptr<IType> iset =
        IType::Instance(csh_s, csh_f_dx, csh_s, csh_s, 0u);
    std::cout << "Created integral set " << iset->label()
              << " key = " << iset->key() << std::endl;
  }
  {
    typedef TwoPRep_11_11_sq IType;
    std::shared_ptr<IType> iset =
        IType::Instance(csh_q, csh_q, csh_q, csh_q, 0u);
    std::cout << "Created integral set " << iset->label()
              << " key = " << iset->key() << std::endl;
  }
}

void test3() {
  {
    typedef TwoPRep_11_11_sq IType;
    typedef VRR_a_11_TwoPRep_11_sh RRType;
    std::shared_ptr<IType> iset = IType::Instance(sh_p, sh_p, sh_p, sh_p, 0u);
    std::shared_ptr<RRType> rr = RRType::Instance(iset, 0);
    std::cout << "Created recurrence relation " << rr->label() << std::endl;
  }
  {
    typedef TwoPRep_11_11_sq IType;
    typedef VRR_c_11_TwoPRep_11_sh RRType;
    std::shared_ptr<IType> iset = IType::Instance(sh_p, sh_p, sh_p, sh_p, 0u);
    std::shared_ptr<RRType> rr = RRType::Instance(iset, 0);
    std::cout << "Created recurrence relation " << rr->label() << std::endl;
  }
  {
    typedef DivG12prime_xTx_11_11_sq IType;
    typedef CR_11_DivG12prime_xTx_11_sh RRType;
    std::shared_ptr<IType> iset =
        IType::Instance(sh_p, sh_p, sh_p, sh_p, 0u, DivG12prime_xTx_Descr(0));
    std::shared_ptr<RRType> rr = RRType::Instance(iset, 0);
    std::cout << "Created recurrence relation " << rr->label() << std::endl;
  }
}

// primitive ERI build
void test4() {
  RunBuildTest<TwoPRep_11_11_sq>(sh_p, sh_s, sh_p, sh_s, 0, use_quartets);
  RunBuildTest<TwoPRep_11_11_sq>(sh_p, sh_p, sh_p, sh_p, 0, use_quartets);
  RunBuildTest<TwoPRep_11_11_sq>(sh_p, sh_p, sh_p, sh_p, 0, use_integrals);
}

// contracted ERI build
void test5() {
  CGShell::set_contracted_default_value(true);
  const bool contracted_targets_old_value = cparams->contracted_targets();
  cparams->contracted_targets(true);
  CGShell csh_s(0u);
  CGShell csh_p(1u);

  // RunBuildTest<TwoPRep_11_11_sq>(csh_p,csh_s,csh_p,csh_s,0,use_quartets);
  RunBuildTest<TwoPRep_11_11_sq>(csh_s, csh_p, csh_s, csh_s, 0, use_quartets);

  cparams->contracted_targets(contracted_targets_old_value);
}

// contracted derivative ERI build
void test6() {
  CGShell::set_contracted_default_value(true);
  const bool contracted_targets_old_value = cparams->contracted_targets();
  cparams->contracted_targets(true);
  CGShell csh_s(0u);
  CGShell csh_p(1u);
  CGShell csh_d(2u);
  CGShell csh_s_dx(0u);
  csh_s_dx.deriv().inc(0, 1);
  CGShell csh_p_dx(1u);
  csh_p_dx.deriv().inc(0, 1);
  CGShell csh_p_dy(1u);
  csh_p_dy.deriv().inc(1, 1);
  CGShell csh_p_dxyz(1u);
  csh_p_dxyz.deriv().inc(0, 1);
  csh_p_dxyz.deriv().inc(1, 1);
  csh_p_dxyz.deriv().inc(2, 1);
  CGShell csh_d_dx(2u);
  csh_d_dx.deriv().inc(0, 1);
  CGShell csh_s_d2x(0u);
  csh_s_d2x.deriv().inc(0, 2);
  CGShell csh_p_d2x(1u);
  csh_p_d2x.deriv().inc(0, 2);

  RunBuildTest<TwoPRep_11_11_sq>(csh_d_dx, csh_s, csh_s, csh_s, 0,
                                 use_integrals);
  // RunBuildTest<TwoPRep_11_11_sq>(csh_d_dx,csh_s,csh_s,csh_s,0,use_quartets);

  cparams->contracted_targets(contracted_targets_old_value);
}

// testing RR generators
void test7() {
  CGShell::set_contracted_default_value(true);
  const bool contracted_targets_old_value = cparams->contracted_targets();
  cparams->contracted_targets(true);

  CGShell csh_s(0u);
  CGShell csh_p(1u);
  CGShell csh_d(2u);
  CGShell csh_s_dx(0u);
  csh_s_dx.deriv().inc(0, 1);
  CGShell csh_p_dx(1u);
  csh_p_dx.deriv().inc(0, 1);
  CGShell csh_p_dy(1u);
  csh_p_dy.deriv().inc(1, 1);
  CGShell csh_p_dxyz(1u);
  csh_p_dxyz.deriv().inc(0, 1);
  csh_p_dxyz.deriv().inc(1, 1);
  csh_p_dxyz.deriv().inc(2, 1);
  CGShell csh_d_dx(2u);
  csh_d_dx.deriv().inc(0, 1);
  CGShell csh_s_d2x(0u);
  csh_s_d2x.deriv().inc(0, 2);
  CGShell csh_p_d2x(1u);
  csh_p_d2x.deriv().inc(0, 2);

  std::shared_ptr<TwoPRep_11_11_sq> target =
      TwoPRep_11_11_sq::Instance(csh_p, csh_d_dx, csh_s, csh_s, mType(0u));

  std::shared_ptr<RecurrenceRelation> rr =
      HRR_ab_11_TwoPRep_11_sh::Instance(target, 0);
  assert(rr != 0 && rr->num_children() != 0);
  std::shared_ptr<RRStack> rrstack = RRStack::Instance();
  rrstack->find(rr);

  std::deque<std::string> decl_filenames, def_filenames;
  generate_rr_code(std::cout, cparams, decl_filenames, def_filenames);

  cparams->contracted_targets(contracted_targets_old_value);
}

// contracted g12 integral build
void test8() {
  CGShell::set_contracted_default_value(true);
  const bool contracted_targets_old_value = cparams->contracted_targets();
  cparams->contracted_targets(true);
  CGShell csh_s(0u);
  CGShell csh_p(1u);

  const Ti_G12_Descr t0g12_descr(0);
  RunBuildTest<TiG12_11_11_sq>(csh_p, csh_s, csh_p, csh_s, 0, t0g12_descr,
                               use_quartets);

  cparams->contracted_targets(contracted_targets_old_value);
}

// DerivGaussV2: the children of direct differentiation and of translational
// invariance, and the labels by which the RRStack shares the generated code
void test9() {
  typedef TwoPRep_11_11_sq IType;
  auto check = [](bool ok, const std::string& what) {
    if (!ok) throw ProgrammingError("DerivGaussV2 test: " + what);
  };
  // @return whether child @p i of @p rr is @p iset
  auto child_is = [](const std::shared_ptr<RecurrenceRelation>& rr,
                     unsigned int i, const std::shared_ptr<IType>& iset) {
    return rr->rr_child(i) == std::static_pointer_cast<DGVertex, IType>(iset);
  };
  // @return shell of angular momentum @p l (primitive, unless @p contracted)
  // with derivative quanta @p dx, @p dy, @p dz
  auto shell = [](unsigned int l, unsigned int dx = 0, unsigned int dy = 0,
                  unsigned int dz = 0, bool contracted = false) {
    CGShell sh(l);
    if (contracted)
      sh.contract();
    else
      sh.uncontract();
    if (dx) sh.deriv().inc(0, dx);
    if (dy) sh.deriv().inc(1, dy);
    if (dz) sh.deriv().inc(2, dz);
    return sh;
  };
  const mType m0(0u);

  // direct differentiation of a (particle 0, bra), in x:
  //   ∂_Ax (p_x s|s s) = 2α (d s|s s) - 1 (s s|s s)
  {
    auto target =
        IType::Instance(shell(1, 1), shell(0), shell(0), shell(0), m0);
    auto rr = DerivV2_a_11_TwoPRep_11_sh::Instance(target, 0);
    check(rr != nullptr, "no relation for (p_x s|s s)");
    check(rr->num_children() == 2, "(p_x s|s s) must have 2 children");
    check(child_is(rr, 0,
                   IType::Instance(shell(2), shell(0), shell(0), shell(0), m0)),
          "child 0 of (p_x s|s s) must be (d s|s s)");
    check(child_is(rr, 1,
                   IType::Instance(shell(0), shell(0), shell(0), shell(0), m0)),
          "child 1 of (p_x s|s s) must be (s s|s s)");
    // s has no a-1 term
    auto target_s =
        IType::Instance(shell(0, 1), shell(0), shell(0), shell(0), m0);
    auto rr_s = DerivV2_a_11_TwoPRep_11_sh::Instance(target_s, 0);
    check(rr_s != nullptr && rr_s->num_children() == 1,
          "(s_x s|s s) must have 1 child");
    check(child_is(rr_s, 0,
                   IType::Instance(shell(1), shell(0), shell(0), shell(0), m0)),
          "child 0 of (s_x s|s s) must be (p s|s s)");
    // no derivative quanta in y, no relation
    check(DerivV2_a_11_TwoPRep_11_sh::Instance(target, 1) == nullptr,
          "(p_x s|s s) must not be differentiated in y");
    // direct differentiation applies to primitives only
    auto target_c =
        IType::Instance(shell(1, 1, 0, 0, true), shell(0, 0, 0, 0, true),
                        shell(0, 0, 0, 0, true), shell(0, 0, 0, 0, true), m0);
    check(DerivV2_a_11_TwoPRep_11_sh::Instance(target_c, 0) == nullptr,
          "contracted (p_x s|s s) must not be differentiated directly");
  }

  // code sharing: the label depends on the differentiated shell and the
  // direction only, not on the spectators
  {
    auto rr1 = DerivV2_a_11_TwoPRep_11_sh::Instance(
        IType::Instance(shell(1, 1), shell(0), shell(0), shell(0), m0), 0);
    auto rr2 = DerivV2_a_11_TwoPRep_11_sh::Instance(
        IType::Instance(shell(1, 1), shell(2), shell(1), shell(1), m0), 0);
    auto rr3 = DerivV2_a_11_TwoPRep_11_sh::Instance(
        IType::Instance(shell(2, 1), shell(0), shell(0), shell(0), m0), 0);
    check(rr1 != nullptr && rr2 != nullptr && rr3 != nullptr,
          "no relation for (p_x s|s s), (p_x d|p p) or (d_x s|s s)");
    check(rr1->label() == rr2->label(),
          "(p_x s|s s) and (p_x d|p p) must share the code: " + rr1->label() +
              " vs " + rr2->label());
    check(rr1->label() != rr3->label(),
          "(p_x s|s s) and (d_x s|s s) must not share the code");
    // a shell with derivative quanta in x and y is differentiated in either
    // direction, with different children: the labels must differ
    auto target_xy =
        IType::Instance(shell(1, 1, 1), shell(0), shell(0), shell(0), m0);
    auto rr_x = DerivV2_a_11_TwoPRep_11_sh::Instance(target_xy, 0);
    auto rr_y = DerivV2_a_11_TwoPRep_11_sh::Instance(target_xy, 1);
    check(rr_x != nullptr && rr_y != nullptr,
          "(p_xy s|s s) must be differentiated in x and y");
    check(rr_x->label() != rr_y->label(),
          "the x and y derivatives of (p_xy s|s s) must not share the code");
  }

  // translational invariance (particle 1, bra = c), also for contracted
  // shells: ∂_Cx (a b|c d) = -∂_Ax (a b|c d) - ∂_Bx (a b|c d) - ∂_Dx (a b|c d)
  {
    typedef DerivGaussV2<IType, 1, InBra, 1, InBra> TIType;
    const bool c = true;
    auto sh = [&](unsigned int l, unsigned int dx = 0) {
      return shell(l, dx, 0, 0, c);
    };
    auto target = IType::Instance(sh(1), sh(0), sh(1, 1), sh(1), m0);
    auto rr = TIType::Instance(target, 0);
    check(rr != nullptr, "no translational invariance relation");
    check(rr->num_children() == 3,
          "translational invariance must have 3 children");
    check(child_is(rr, 0, IType::Instance(sh(1, 1), sh(0), sh(1), sh(1), m0)),
          "child 0 must be (p_x s|p p)");
    check(child_is(rr, 1, IType::Instance(sh(1), sh(0, 1), sh(1), sh(1), m0)),
          "child 1 must be (p s_x|p p)");
    check(child_is(rr, 2, IType::Instance(sh(1), sh(0), sh(1), sh(1, 1), m0)),
          "child 2 must be (p s|p p_x)");
    check(rr->label().rfind("DerivGaussV2TI", 0) == 0,
          "translational invariance label must start with DerivGaussV2TI: " +
              rr->label());
  }

  // build contracted shell sets, without unrolling, hence with calls to the
  // generated DerivGaussV2 functions: differentiation of a (direct, shared by
  // the spectators) and of c (translational invariance)
  {
    CGShell::set_contracted_default_value(true);
    const bool contracted_targets_old_value = cparams->contracted_targets();
    cparams->contracted_targets(true);
    const bool c = true;
    RunBuildTest<IType>(shell(1, 1, 0, 0, c), shell(1, 0, 0, 0, c),
                        shell(0, 0, 0, 0, c), shell(1, 0, 0, 0, c), 0,
                        use_quartets);
    RunBuildTest<IType>(shell(1, 0, 0, 0, c), shell(0, 0, 0, 0, c),
                        shell(1, 0, 1, 0, c), shell(0, 0, 0, 0, c), 0,
                        use_quartets);
    cparams->contracted_targets(contracted_targets_old_value);
  }

  // CR_11_Coulombσp_11 is an alias (target = child, see
  // DirectedGraph::handle_trivial_nodes()), CR_11_Coulombσpσp_11 is not
  {
    for (int k = 0; k != 3; ++k) {
      auto target = Coulombσp_11_11_sq::Instance(
          shell(1), shell(0), shell(1), shell(1), m0, Coulombσp_Descr(k));
      auto rr = CR_11_Coulombσp_11_sh::Instance(target, 0);
      check(rr != nullptr && rr->num_children() == 1 && rr->is_alias(),
            "CR_11_Coulombσp_11 must be an alias");
    }
  }
}

};  // namespace
