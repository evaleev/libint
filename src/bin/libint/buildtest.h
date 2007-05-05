
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

  /// Command-line parser for the standard build tester -- N is the number of centers, i.e. 4 for 4-center ERI
  template <unsigned int N>
    class TesterCmdLine{
  public:
    TesterCmdLine(int argc, char* argv[]);
    ~TesterCmdLine() {}
    
    const std::vector<unsigned int>& am() const { return am_; }
    unsigned int size_to_unroll() const { return size_to_unroll_; }
    unsigned int veclen() const { return veclen_; }
    bool vectorize_by_line() const { return vectorize_by_line_; }
    bool do_cse() const { return do_cse_; }

  private:
    static const unsigned int max_am = 10;
    std::vector<unsigned int> am_;
    unsigned int size_to_unroll_;
    unsigned int veclen_;
    bool vectorize_by_line_;
    bool do_cse_;
  };

  /** This is a user-friendly generic test of building an Integral using specified size_to_unroll, veclen, vec_by_line, and do_cse.
      GenAllCode should be set to true if compilable code
      to be produced (i.e. include header files + set-level recurrence relations code)
   */
  template <class Integral, bool GenAllCode>
    void BuildTest(const SafePtr<Integral>& target, unsigned int size_to_unroll, unsigned int veclen,
		   bool vec_by_line, bool do_cse, const std::string& complabel = "buildtest",
		   std::ostream& os = std::cout);

  /** This is a generic test of building an Integral using specified cparams, memman, size_to_unroll,
      default strategy and specified tactic. GenAllCode should be set to true if compilable code
      to be produced (i.e. include header files + set-level recurrence relations code)
   */
  template <class Integral, bool GenAllCode>
    void __BuildTest(const SafePtr<Integral>& target, const SafePtr<CompilationParameters>& cparams,
		     unsigned int size_to_unroll, std::ostream& os = std::cout,
		     const SafePtr<Tactic>& tactic = SafePtr<Tactic>(new FirstChoiceTactic<DummyRandomizePolicy>),
		     const SafePtr<MemoryManager>& memman = SafePtr<MemoryManager>(new WorstFitMemoryManager),
		     const std::string& complabel = "general_integral");

  template <class Integral, bool GenAllCode>
    void
    __BuildTest(const SafePtr<Integral>& target, const SafePtr<CompilationParameters>& cparams,
		unsigned int size_to_unroll, std::ostream& os,
		const SafePtr<Tactic>& tactic, const SafePtr<MemoryManager>& memman,
		const std::string& complabel)
    {
      const std::string label = cparams->api_prefix() + target->label();
      SafePtr<DirectedGraph> dg_xxxx(new DirectedGraph);
      SafePtr<Strategy> strat(new Strategy(size_to_unroll));
      os << "Building " << target->description() << std::endl;
      SafePtr<DGVertex> xsxs_ptr = dynamic_pointer_cast<DGVertex,Integral>(target);

      LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
      taskmgr.add(complabel);
      taskmgr.current(complabel);
      
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
      const bool need_to_optimize = (max_am <= cparams->max_am_opt(complabel));
      dg_xxxx->registry()->do_cse(need_to_optimize);
      dg_xxxx->registry()->condense_expr(condense_expr(cparams->unroll_threshold(),cparams->max_vector_length()>1));
      
      // Need to accumulate integrals?
      dg_xxxx->registry()->accumulate_targets(cparams->accumulate_targets());

      // Only build using shell-sets?
      if (size_to_unroll == 0)
	dg_xxxx->registry()->can_unroll(false);

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
      taskmgr.current().params()->max_stack_size(memman->max_memory_used());
      // extract all extrnal symbols
      extract_symbols(dg_xxxx);

      std::basic_ofstream<char> dotfile2("graph.symb.dot");
      dg_xxxx->print_to_dot(true,dotfile2);
      
      os << "Max memory used = " << memman->max_memory_used() << std::endl;
      dg_xxxx->reset();

      if (GenAllCode) {
	// initialize code context to produce library API
	SafePtr<CodeContext> icontext(new CppCodeContext(cparams));
	// initialize object to generate interface
	SafePtr<Libint2Iface> iface(new Libint2Iface(cparams,icontext));

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


  template <class Integral, bool GenAllCode>
    void BuildTest(const SafePtr<Integral>& target, unsigned int size_to_unroll, unsigned int veclen,
		   bool vec_by_line, bool do_cse, const std::string& complabel,
		   std::ostream& os)
  {
    const unsigned int max_am = 10;
    os << "generating code to compute " << target->label() << std::endl;

    LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
    taskmgr.add(complabel);
    taskmgr.current(complabel);
      
    // initialize cparams
    SafePtr<CompilationParameters> cparams(new CompilationParameters);
    cparams->max_am(complabel,max_am);
    cparams->max_vector_length(veclen);
    cparams->vectorize_by_line(vec_by_line);
    cparams->count_flops(true);
#if LIBINT_ACCUM_INTS
    cparams->accumulate_targets(true);
#else
    cparams->accumulate_targets(false);
#endif
#ifdef LIBINT_API_PREFIX
  {
    const std::string api_prefix(LIBINT_API_PREFIX);
    cparams->api_prefix(api_prefix);
  }
#endif

    if (do_cse) {
      cparams->max_am_opt(complabel,max_am);
    }
    else {
      cparams->max_am_opt(complabel,0);
    }

    // set default dims
    ImplicitDimensions::set_default_dims(cparams);

    SafePtr<StdRandomizePolicy> rpolicy(new StdRandomizePolicy(0.00));
    SafePtr<Tactic> tactic(new FirstChoiceTactic<StdRandomizePolicy>(rpolicy));
    const SafePtr<MemoryManager> memman(new WorstFitMemoryManager);
    __BuildTest<Integral,true>(target,cparams,size_to_unroll,os,tactic,memman,complabel);
  }

  template <unsigned int N>
    TesterCmdLine<N>::TesterCmdLine(int argc, char* argv[])
    {
      if (N == 0)
	throw ProgrammingError("TesterCmdLine<N>::TesterCmdLine but N is 0");
      const int argc_min = N + 2;
      const int argc_max = N + 5;
      if (argc < argc_min || argc > argc_max) {
	std::cerr << "Usage: " << argv[0] << " <am> size_to_unroll [vector_length] [vector_method] [do_cse]" << std::endl
		  << "       <am> -- angular momenta on each center, e.g. 4 nonnegative integers for a 4-center ERI" << std::endl
		  << "       size_to_unroll -- size of the largest integral set to be unrolled" << std::endl
		  << "       vector_length  -- (optional) max vector length. Defaults to 1." << std::endl
		  << "       vector_method  -- (optional) vectorization method. Valid choices are 0 (by-block) and 1 (by-line). Defaults to 0." << std::endl
		  << "       do_cse  -- (optional) do Common Subexpression Elimination? Valid choices are 0 (no) and 1 (yes). Defaults to 0." << std::endl << std::endl;
	throw InputError("TesterCmdLine<N>::TesterCmdLine -- incorrect number of command-line arguments");
      }
      for(unsigned int i=1; i<N+1; ++i) {
	const unsigned int am = atoi(argv[i]);
	if (am > max_am)
	  throw InputError("TesterCmdLine<N>::TesterCmdLine -- angular momentum limit exceeded");
	am_.push_back(am);
      }
      size_to_unroll_ = atoi(argv[N+1]);

      veclen_ = 1;
      if (argc >= N+3) {
	veclen_ = atoi(argv[N+2]);
      }
      vectorize_by_line_ = false;
      if (argc >= N+4) {
	vectorize_by_line_ = (1 == atoi(argv[N+3]));
      }
      do_cse_ = false;
      if (argc >= N+5) {
	do_cse_ = (1 == atoi(argv[N+4]));
      }
    }

};

#endif
