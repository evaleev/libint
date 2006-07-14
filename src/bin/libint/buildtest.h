
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <dg.h>
#include <strategy.h>
#include <iface.h>
#include <graph_registry.h>

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
		   const SafePtr<Tactic>& tactic = SafePtr<Tactic>(new FirstChoiceTactic<DummyRandomizePolicy>),
		   const SafePtr<MemoryManager>& memman = SafePtr<MemoryManager>(new WorstFitMemoryManager),
		   const std::string& complabel = "general_integral");

  template <class Integral, bool GenAllCode>
    void
    BuildTest(const SafePtr<Integral>& target, const SafePtr<CompilationParameters>& cparams,
	      unsigned int size_to_unroll, std::ostream& os,
	      const SafePtr<Tactic>& tactic, const SafePtr<MemoryManager>& memman,
	      const std::string& complabel)
    {
      const std::string label = cparams->api_prefix() + target->label();
      SafePtr<DirectedGraph> dg_xxxx(new DirectedGraph);
      SafePtr<Strategy> strat(new Strategy(size_to_unroll));
      os << "Building " << target->description() << std::endl;
      SafePtr<DGVertex> xsxs_ptr = dynamic_pointer_cast<DGVertex,Integral>(target);
      
      //
      // do CSE only if max_am <= cparams->max_am_opt()
      //
      const unsigned int np = target->bra().num_part();
      unsigned int max_am = 0;
      // bra
      for(unsigned int p=0; p<np; p++) {
        const unsigned int nf = target->bra().num_members(p);
        for(unsigned int f=0; f<nf; f++) {
          // Assuming shells here
          const unsigned int am = target->bra(p,f).qn();
          using std::max;
          max_am = max(max_am,am);
        }
      }
      // ket
      for(unsigned int p=0; p<np; p++) {
        const unsigned int nf = target->ket().num_members(p);
        for(unsigned int f=0; f<nf; f++) {
          // Assuming shells here
          const unsigned int am = target->ket(p,f).qn();
          using std::max;
          max_am = max(max_am,am);
        }
      }
      const bool need_to_optimize = (max_am <= cparams->max_am_opt());
      dg_xxxx->registry()->do_cse(need_to_optimize);
      
      // Need to accumulate integrals?
      dg_xxxx->registry()->accumulate_targets(cparams->accumulate_targets());

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
      
      // update max stack size
      LibraryParameters& lparams = LibraryParameters::get_library_params();
      lparams.max_stack_size(memman->max_memory_used());
      
      std::basic_ofstream<char> dotfile2("graph.symb.dot");
      dg_xxxx->print_to_dot(true,dotfile2);
      
      os << "Max memory used = " << memman->max_memory_used() << std::endl;
      dg_xxxx->reset();

      if (GenAllCode) {
	// initialize code context to produce library API
	SafePtr<CodeContext> icontext(new CppCodeContext(cparams));
	// make a list of computation labels
	Libint2Iface::Comps comps;
	comps.push_back(complabel);
	// initialize object to generate interface
	SafePtr<Libint2Iface> iface(new Libint2Iface(cparams,icontext,comps));

	// generate interface
	std::ostringstream oss;
	oss << "#include <" << decl_filename << ">" << std::endl;
	iface->to_int_iface(oss.str());

	// Generate set-level RR code
	generate_rr_code(os,cparams);

	// Print log
	std::cout << "Generated headers: " << decl_filename;
	SafePtr<RRStack> rrstack = RRStack::Instance();
	for(RRStack::citer_type it = rrstack->begin(); it!=rrstack->end(); it++) {
	  SafePtr<RecurrenceRelation> rr = (*it).second.second;
	  std::string rrlabel = cparams->api_prefix() + rr->label();
	  std::cout << " " << context->label_to_name(rrlabel) << ".h";
	}
	std::cout << std::endl << "Generated sources: " << def_filename;
	for(RRStack::citer_type it = rrstack->begin(); it!=rrstack->end(); it++) {
	  SafePtr<RecurrenceRelation> rr = (*it).second.second;
	  std::string rrlabel = cparams->api_prefix() + rr->label();
	  std::cout << " " << context->label_to_name(rrlabel) << ".cc";
	}
	std::cout << std::endl << "Top compute function: " << context->label_to_name(label_to_funcname(label)) << std::endl;

      }
    }


};

#endif
