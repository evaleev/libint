
#include <string>
#include <fstream>
#include <iface.h>
#include <default_params.h>

using namespace std;
using namespace libint2;

namespace {
  const char mh_name[] = "libint2.h";
  const char th_name[] = "libint2_types.h";
  const char ph_name[] = "libint2_params.h";
  const char ih_name[] = "libint2_iface.h";
  const char ii_name[] = "libint2_iface_internal.h";
  const char si_name[] = "libint2_static_init.cc";
  const char sc_name[] = "libint2_static_cleanup.cc";
  const char li_name[] = "libint2_init.cc";
  
  inline void header_guard_open(std::ostream& os, const std::string& label) {
    os << "#ifndef _libint2_" << label << "_h_" << endl
       << "#define _libint2_" << label << "_h_" << endl << endl;
  }
  inline void header_guard_close(std::ostream& os) {
    os << "#endif" << endl << endl;
  }

};

Libint2Iface::Libint2Iface(const SafePtr<CompilationParameters>& cparams,
                           const SafePtr<CodeContext>& ctext) :
  cparams_(cparams), ctext_(ctext), null_str_(""),
  th_((cparams_->source_directory() + th_name).c_str()),
  ph_((cparams_->source_directory() + ph_name).c_str()),
  ih_((cparams_->source_directory() + ih_name).c_str()),
  ii_((cparams_->source_directory() + ii_name).c_str()),
  si_((cparams_->source_directory() + si_name).c_str()),
  sc_((cparams_->source_directory() + sc_name).c_str()),
  li_((cparams_->source_directory() + li_name).c_str())
{
  header_guard_open(th_,"libint2types");
  header_guard_open(ph_,"libint2params");
  header_guard_open(ih_,"libint2iface");
  header_guard_open(ii_,"libint2ifaceint");

  ph_ << macro_define("API_PREFIX", cparams_->api_prefix());
  ph_ << macro_define("MAX_VECLEN",cparams_->max_vector_length());
  if (cparams_->count_flops())
    ph_ << macro_define("FLOP_COUNT",1);
  if (cparams_->accumulate_targets())
    ph_ << macro_define("ACCUM_INTS",1);
  const std::string realtype(cparams_->realtype());
  ph_ << macro_define("REALTYPE",realtype);
  
  ih_ << "#include <cstddef>" << endl
      << ctext_->code_prefix();

  oss_.str(null_str_);
  oss_ << ctext_->std_header() << "#include <" << ih_name << ">" << endl
                               << "#include <" << ii_name << ">" << endl
                               << "#include <cstddef>" << endl
                               << ctext_->code_prefix();
  std::string pfix = oss_.str();
  si_ << pfix;
  sc_ << pfix;
  li_ << pfix;

  // print out declarations for the array of pointers to evaluator functions
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  typedef LibraryTaskManager::TasksCIter tciter;
  for(tciter t=taskmgr.first(); t!=taskmgr.plast(); ++t) {
    const std::string& tlabel = t->label();
    const unsigned int lmax = cparams_->max_am(tlabel) + 1;

    ostringstream oss;
    oss << "void (*" << ctext->label_to_name(cparams->api_prefix()) << "libint2_build_" << tlabel
        << "[" << lmax << "][" << lmax << "][" << lmax << "]["
        << lmax << "])(" << ctext_->const_modifier() << ctext_->inteval_type_name(tlabel) << "*);" << endl;
    ih_ << "extern " << oss.str();
    si_ << oss.str();
  }
  
  // Declare library constructor/destructor
  oss_.str(null_str_);
  oss_ << ctext_->type_name<void>() << " "
       << ctext_->label_to_name(cparams->api_prefix() + "libint2_static_init") << "()";
  std::string si_fdec(oss_.str());
  
  oss_.str(null_str_);
  oss_ << ctext_->type_name<void>() << " "
       << ctext_->label_to_name(cparams->api_prefix() + "libint2_static_cleanup") << "()";
  std::string sc_fdec(oss_.str());

  ih_ << si_fdec << ctext_->end_of_stat() << endl;
  ih_ << sc_fdec << ctext_->end_of_stat() << endl;

  // Declare evaluator constructor/destructor (specific to the type of computation)
  for(tciter t=taskmgr.first(); t!=taskmgr.plast(); ++t) {
    const std::string& tlabel = t->label();

    oss_.str(null_str_);
    oss_ << ctext_->type_name<void>() << " "
	 << ctext_->label_to_name(cparams->api_prefix() + "libint2_init_" + tlabel) << "(" << ctext_->inteval_type_name(tlabel) << "* inteval, int max_am, LIBINT2_REALTYPE* buf)";
    std::string li_fdec(oss_.str());
    li_decls_.push_back(li_fdec);
  
    oss_.str(null_str_);
    oss_ << ctext_->type_name<size_t>() << " "
	 << ctext_->label_to_name(cparams->api_prefix() + "libint2_need_memory_" + tlabel) << "(int max_am)";
    std::string lm_fdec(oss_.str());
    lm_decls_.push_back(lm_fdec);

    oss_.str(null_str_);
    oss_ << ctext_->type_name<void>() << " "
	 << ctext_->label_to_name(cparams->api_prefix() + "libint2_cleanup_" + tlabel) << "(" << ctext_->inteval_type_name(tlabel) << "* inteval)";
    std::string lc_fdec(oss_.str());
    lc_decls_.push_back(lc_fdec);

    ih_ << li_fdec << ctext_->end_of_stat() << endl;
    ih_ << lm_fdec << ctext_->end_of_stat() << endl;
    ih_ << lc_fdec << ctext_->end_of_stat() << endl;
  }

  ih_ << ctext_->code_postfix() << endl;
  
  si_ << si_fdec << ctext_->open_block();
  sc_ << sc_fdec << ctext_->open_block();
  
}

Libint2Iface::~Libint2Iface()
{
  // For each task, print out defines for stack dimensions
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  typedef LibraryTaskManager::TasksCIter tciter;
  for(tciter t=taskmgr.first(); t!=taskmgr.plast(); ++t) {
    SafePtr<TaskParameters> tparams = t->params();
    const std::string& tlabel = t->label();
    ph_ << macro_define(tlabel,"NUM_TARGETS",tparams->max_ntarget());
    ph_ << macro_define(tlabel,"MAX_STACK_SIZE",tparams->max_stack_size());
    const unsigned int max_vector_stack_size = tparams->max_vector_stack_size();
    if (max_vector_stack_size)
      ph_ << macro_define(tlabel,"MAX_VECTOR_STACK_SIZE",max_vector_stack_size);
    ph_ << macro_define(tlabel,"MAX_HRR_HSRANK",tparams->max_hrr_hsrank());
    ph_ << macro_define(tlabel,"MAX_HRR_LSRANK",tparams->max_hrr_lsrank());
  }

  // For each task, generate the evaluator type
  generate_inteval_type(th_);

  // libint2_iface.h needs macros to help forming prefixed names in API
  ih_ << ctext_->comment("Use LIBINT2_PREFIXED_NAME(fncname) to form properly prefixed function name from LIBINT2 API") << std::endl;
  ih_ << "#define LIBINT2_PREFIXED_NAME(name) __libint2_prefixed_name__(LIBINT2_API_PREFIX,name)" << std::endl;
  ih_ << "#define __libint2_prefixed_name__(prefix,name) __prescanned_prefixed_name__(prefix,name)" << std::endl;
  ih_ << "#define __prescanned_prefixed_name__(prefix,name) prefix##name" << std::endl;
  // also need macros to test what's in the evaluator
  ih_ << ctext_->comment("Use LIBINT2_PREFIXED_NAME(fncname) to form properly prefixed function name from LIBINT2 API") << std::endl;
  ih_ << "#define LIBINT2_DEFINED(taskname,symbol) __prescanned_libint2_defined__(taskname,symbol)" << std::endl;
  if (cparams_->single_evaltype())
    ih_ << "#define __prescanned_libint2_defined__(taskname,symbol) LIBINT2_DEFINED_##symbol" << std::endl << std::endl;
  else
    ih_ << "#define __prescanned_libint2_defined__(taskname,symbol) LIBINT2_DEFINED_##symbol##_##taskname" << std::endl << std::endl;


  
  header_guard_close(th_);
  header_guard_close(ph_);
  header_guard_close(ih_);
  header_guard_close(ii_);

  std::ostringstream oss;
  oss << ctext_->close_block() << ctext_->code_postfix() << endl << endl;
  std::string pfix = oss.str();

  si_ << pfix;
  sc_ << pfix;

  // Define Libint_t constructor/destructor (specific to the type of computation)
  {
    unsigned int i = 0;
    for(tciter t=taskmgr.first(); t!=taskmgr.plast(); ++t,++i) {
      const std::string& tlabel = t->label();

      li_ << lm_decls_[i] << ctext_->open_block();
      li_ << "return " << macro(tlabel,"MAX_STACK_SIZE") << " * " << macro("MAX_VECLEN") << " + "
	  << macro(tlabel,"MAX_VECTOR_STACK_SIZE") << " * " << macro("MAX_VECLEN") << " * ("
	    << macro(tlabel,"MAX_HRR_HSRANK") << " > " << macro(tlabel,"MAX_HRR_LSRANK") << " ? "
	      << macro(tlabel,"MAX_HRR_HSRANK") << " : " << macro(tlabel,"MAX_HRR_LSRANK") << ");" << std::endl;
      li_ << ctext_->close_block();
    }

    i = 0;
    for(tciter t=taskmgr.first(); t!=taskmgr.plast(); ++t,++i) {
      const std::string& tlabel = t->label();

      li_ << li_decls_[i] << ctext_->open_block();
      li_ << "if (buf != 0) inteval->stack = buf;" << std::endl << "else ";
      {
	std::string tmp = ctext_->label_to_name(cparams_->api_prefix() + "libint2_need_memory_" + tlabel) + "(max_am)";
	li_ << ctext_->assign("inteval->stack","new LIBINT2_REALTYPE[" + tmp + "]");
      }

      std::string vstack_ptr("inteval->stack + ");
      vstack_ptr += macro(tlabel,"MAX_STACK_SIZE");
      vstack_ptr += " * ";
      vstack_ptr += macro("MAX_VECLEN");

      li_ << ctext_->assign("inteval->vstack",vstack_ptr);
      if (cparams_->count_flops()) {
	// set the counter to zero
	li_ << "inteval->nflops = 0;" << endl;
      }
      li_ << ctext_->close_block();
    }

    i = 0;
    for(tciter t=taskmgr.first(); t!=taskmgr.plast(); ++t,++i) {
      li_ << lc_decls_[i] << ctext_->open_block();
      li_ << "delete[] inteval->stack;" << std::endl;
      li_ << ctext_->assign("inteval->stack","0");
      li_ << ctext_->assign("inteval->vstack","0");
      if (cparams_->count_flops()) {
	// set the counter to zero
	li_ << "inteval->nflops = 0;" << endl;
      }
      li_ << ctext_->close_block();
    }
  }

  li_ << ctext_->code_postfix() << std::endl << std::endl;

  th_.close();
  ph_.close();
  ih_.close();
  ii_.close();
  si_.close();
  sc_.close();
  li_.close();
}


void
Libint2Iface::generate_inteval_type(std::ostream& os)
{
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();

  //
  // If need to generate single type for all tasks, take a union of all symbols
  // else process each task separately
  //
  typedef LibraryTaskManager::TasksCIter tciter;
  const tciter tend = cparams_->single_evaltype() ? taskmgr.first()+1 : taskmgr.plast();
  for(tciter t=taskmgr.first(); t!=tend; ++t) {
    const SafePtr<TaskExternSymbols> tsymbols = t->symbols();
    
    // Prologue
    os << "typedef struct {" << std::endl;

    //
    // Declare external symbols
    //
    typedef TaskExternSymbols::SymbolList SymbolList;
    std::string tlabel;
    SymbolList symbols;
    if (cparams_->single_evaltype()) {
      TaskExternSymbols composite_symbols;
      const tciter tend = taskmgr.plast();
      for(tciter t=taskmgr.first(); t!=tend; ++t) {
	const SafePtr<TaskExternSymbols> tsymbols = t->symbols();
	composite_symbols.add(tsymbols->symbols());
      }
      symbols = composite_symbols.symbols();
      tlabel = "";
    }
    else {
      symbols = tsymbols->symbols();
      tlabel = t->label();
    }
    typedef SymbolList::const_iterator citer;
    citer end = symbols.end();
    for(citer s=symbols.begin(); s!=end; ++s) {
      // for each extrnal symbol #define a macro
      std::string tmplabel("DEFINED_"); tmplabel += ctext_->label_to_name(*s);
      os << macro_define(tlabel,tmplabel,1);
      // all symbols are doubles
      os << var_declare_v<double>(*s);
    }

    // Declare members common to all evaluators
    os << ctext_->comment("Scratch buffer to hold intermediates") << std::endl;
    os << ctext_->declare(ctext_->mutable_modifier() + ctext_->type_name<double*>(),std::string("stack"));

    os << ctext_->comment("Buffer to hold vector intermediates. Only used by set-level RR code if it is vectorized linewise") << std::endl;
    os << ctext_->declare(ctext_->mutable_modifier() + ctext_->type_name<double*>(),std::string("vstack"));

    os << ctext_->comment("On completion, this contains pointers to computed targets") << std::endl;
    if (cparams_->single_evaltype()) {
      // figure out the maximum number of targets
      unsigned int max_ntargets = 0;
      const tciter tend = taskmgr.plast();
      for(tciter t=taskmgr.first(); t!=tend; ++t) {
	SafePtr<TaskParameters> tparams = t->params();
	max_ntargets = std::max(max_ntargets,tparams->max_ntarget());
      }
      ostringstream oss;
      oss << max_ntargets;
      os << ctext_->declare_v(ctext_->mutable_modifier() + ctext_->type_name<double*>(),std::string("targets"),oss.str());
    }
    else {
      os << ctext_->declare_v(ctext_->mutable_modifier() + ctext_->type_name<double*>(),std::string("targets"),macro(tlabel,"NUM_TARGETS"));
    }

    os << ctext_->comment("Actual vector length. Not to exceed MAX_VECLEN! If MAX_VECLEN is 1 then veclen is not used") << std::endl;
    os << ctext_->declare(ctext_->type_name<int>(),std::string("veclen"));

    os << ctext_->macro_if(macro("FLOP_COUNT"));
    os << ctext_->comment("FLOP counter. Libint must be configured with --enable-flop-counter to allow FLOP counting. It is user's reponsibility to set zero nflops before computing integrals.") << std::endl;
    os << ctext_->declare(ctext_->mutable_modifier() + macro("UINT_LEAST64"),std::string("nflops"));
    os << ctext_->macro_endif();

    os << ctext_->macro_if(macro("ACCUM_INTS"));
    os << ctext_->comment("If libint was configured with --enable-accum-ints then the target integrals are accumulated. To zero out.the targets automatically before the computation, set this to nonzero.") << std::endl;
    os << ctext_->declare(ctext_->type_name<int>(),std::string("zero_out_targets"));
    os << ctext_->macro_endif();
    
    // Epilogue
    os << "} " << ctext_->inteval_type_name(tlabel) << ctext_->end_of_stat() << std::endl;
    
  }

  // If generating single evaluator type, create aliases from the specialized types to the actual type
  if (cparams_->single_evaltype()) {
    const tciter tend = taskmgr.plast();
    for(tciter t=taskmgr.first(); t!=tend; ++t) {
      const std::string& tlabel = t->label();
      os << "typedef "
	 << ctext_->inteval_gen_type_name() << " "
	 << ctext_->inteval_spec_type_name(tlabel)
	 << ctext_->end_of_stat() << std::endl;
    }
  }

}
