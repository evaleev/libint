
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <dg.h>
#include <strategy.h>
#include <iface.h>

#ifndef _libint2_src_bin_libint_buildtest_h_
#define _libint2_src_bin_libint_buildtest_h_

namespace libint2 {

  // defined in buildtest.cc
  void generate_rr_code(std::ostream& os, const SafePtr<CompilationParameters>& cparams);

  /** This is a generic test of building an Integral using specified cparams, memman, size_to_unroll,
      default strategy and specified tactic. GenAllCode should be set to true if compilable code
      to be produced (i.e. include header files + set-level recurrence relations code)
   */
  template <class Integral, bool GenAllCode>
    void BuildTest(const SafePtr<Integral>& target, const SafePtr<CompilationParameters>& cparams,
		   unsigned int size_to_unroll, std::ostream& os = std::cout,
		   const SafePtr<Tactic>& tactic = SafePtr<Tactic>(new FirstChoiceTactic),
		   const SafePtr<MemoryManager>& memman = SafePtr<MemoryManager>(new WorstFitMemoryManager));

  template <class Integral, bool GenAllCode>
    void
    BuildTest(const SafePtr<Integral>& target, const SafePtr<CompilationParameters>& cparams,
	      unsigned int size_to_unroll, std::ostream& os,
	      const SafePtr<Tactic>& tactic, const SafePtr<MemoryManager>& memman)
    {
      const std::string label = target->label();
      SafePtr<DirectedGraph> dg_xxxx(new DirectedGraph);
      SafePtr<Strategy> strat(new Strategy(size_to_unroll));
      os << "Building " << target->description() << std::endl;
      SafePtr<DGVertex> xsxs_ptr = dynamic_pointer_cast<DGVertex,Integral>(target);
      dg_xxxx->append_target(xsxs_ptr);
      dg_xxxx->apply(strat,tactic);
      dg_xxxx->optimize_rr_out();
      
      std::basic_ofstream<char> dotfile("graph.dot");
      dg_xxxx->print_to_dot(false,dotfile);
      os << "The number of vertices = " << dg_xxxx->num_vertices() << std::endl;
      
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
      
      os << "Max memory used = " << memman->max_memory_used() << std::endl;
      dg_xxxx->reset();

      if (GenAllCode) {
	// initialize code context to produce library API
	SafePtr<CodeContext> icontext(new CppCodeContext(cparams));
	// make a list of computation labels
	Libint2Iface::Comps comps;
	comps.push_back("general_integral");
	// initialize object to generate interface
	SafePtr<Libint2Iface> iface(new Libint2Iface(cparams,icontext,comps));

	// generate interface
	std::ostringstream oss;
	oss << "#include <" << decl_filename << ">" << std::endl;
	iface->to_int_iface(oss.str());

	// Generate set-level RR code
	generate_rr_code(os,cparams);

	// Print log
	std::cout << "Generated header: " << decl_filename << std::endl;
	std::cout << "Generated sources: " << def_filename;
	SafePtr<RRStack> rrstack = RRStack::Instance();
	for(RRStack::citer_type it = rrstack->begin(); it!=rrstack->end(); it++) {
	  SafePtr<RecurrenceRelation> rr = (*it).second.second;
	  std::string rrlabel = rr->label();
	  std::cout << " " << context->label_to_name(rrlabel) << ".cc";
	}
	std::cout << std::endl << "Top compute function: compute" << context->label_to_name(label) << std::endl;

      }
    }


};

#endif
