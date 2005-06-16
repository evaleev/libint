
#include <iostream>
#include <fstream>
#include <string>
#include <dg.h>
#include <strategy.h>

#ifndef _libint2_src_bin_libint_buildtest_h_
#define _libint2_src_bin_libint_buildtest_h_

namespace libint2 {

  // This is a generic test of building an Integral using specified cparams, memman, size_to_unroll, default strategy and specified tactic
  template <class Integral>
    void BuildTest(const SafePtr<Integral>& target, const SafePtr<CompilationParameters>& cparams,
		   unsigned int size_to_unroll, std::ostream& os = std::cout,
		   const SafePtr<Tactic>& tactic = SafePtr<Tactic>(new FirstChoiceTactic),
		   const SafePtr<MemoryManager>& memman = SafePtr<MemoryManager>(new WorstFitMemoryManager));

  template <class Integral>
    void
    BuildTest(const SafePtr<Integral>& target, const SafePtr<CompilationParameters>& cparams,
	      unsigned int size_to_unroll, std::ostream& os,
	      const SafePtr<Tactic>& tactic, const SafePtr<MemoryManager>& memman)
    {
      const std::string label = target->label();
      SafePtr<DirectedGraph> dg_xxxx(new DirectedGraph);
      SafePtr<Strategy> strat(new Strategy(size_to_unroll));
      os << "Building " << target->description() << endl;
      SafePtr<DGVertex> xsxs_ptr = dynamic_pointer_cast<DGVertex,Integral>(target);
      dg_xxxx->append_target(xsxs_ptr);
      dg_xxxx->apply(strat,tactic);
      dg_xxxx->optimize_rr_out();
      
      std::basic_ofstream<char> dotfile("graph.dot");
      dg_xxxx->print_to_dot(false,dotfile);
      os << "The number of vertices = " << dg_xxxx->num_vertices() << endl;
      
      dg_xxxx->traverse();
      //dg_xxxx->debug_print_traversal(cout);
      SafePtr<CodeContext> context(new CppCodeContext(cparams));
      std::string decl_filename(context->label_to_name(label));  decl_filename += ".h";
      std::string def_filename(context->label_to_name(label));  def_filename += ".cc";
      std::basic_ofstream<char> declfile(decl_filename.c_str());
      std::basic_ofstream<char> deffile(def_filename.c_str());
      dg_xxxx->generate_code(context,memman,ImplicitDimensions::default_dims(),SafePtr<CodeSymbols>(new CodeSymbols),
			     label,declfile,deffile);
      
      std::basic_ofstream<char> dotfile2("graph.symb.dot");
      dg_xxxx->print_to_dot(true,dotfile2);
      
      os << "Max memory used = " << memman->max_memory_used() << endl;
      dg_xxxx->reset();
    }


};

#endif
