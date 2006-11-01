
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

  // constructs type name for the evaluator given computation label, e.g. eri -> Libint_eri_t
  inline std::string to_eval_type(const std::string& clabel) {
    std::ostringstream oss;
    oss << "Libint_" << clabel << "_t";
    return oss.str();
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

  ph_ << define("API_PREFIX", cparams_->api_prefix());
  ph_ << define("MAX_VECLEN",cparams_->max_vector_length());
  if (cparams_->count_flops())
    ph_ << define("FLOP_COUNT",1);
  if (cparams_->accumulate_targets())
    ph_ << define("ACCUM_INTS",1);
  const std::string realtype(cparams_->realtype());
  ph_ << define("REALTYPE",realtype);
  
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
  for(tciter t=taskmgr.first(); t!=taskmgr.last(); ++t) {
    const std::string& tlabel = t->label();
    const unsigned int lmax = cparams_->max_am_eri() + 1;

    ostringstream oss;
    oss << "void (*" << ctext->label_to_name(cparams->api_prefix()) << "libint2_build_" << tlabel
        << "[" << lmax << "][" << lmax << "][" << lmax << "]["
        << lmax << "])(" << to_eval_type(tlabel) << "*);" << endl;
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
  for(tciter t=taskmgr.first(); t!=taskmgr.last(); ++t) {
    const std::string& tlabel = t->label();

    oss_.str(null_str_);
    oss_ << ctext_->type_name<void>() << " "
	 << ctext_->label_to_name(cparams->api_prefix() + "libint2_init_" + tlabel) << "(" << to_eval_type(tlabel) << "* inteval, int max_am, LIBINT2_REALTYPE* buf)";
    std::string li_fdec(oss_.str());
    li_decls_.push_back(li_fdec);
  
    oss_.str(null_str_);
    oss_ << ctext_->type_name<size_t>() << " "
	 << ctext_->label_to_name(cparams->api_prefix() + "libint2_need_memory_" + tlabel) << "(int max_am)";
    std::string lm_fdec(oss_.str());
    lm_decls_.push_back(lm_fdec);

    oss_.str(null_str_);
    oss_ << ctext_->type_name<void>() << " "
	 << ctext_->label_to_name(cparams->api_prefix() + "libint2_cleanup_" + tlabel) << "(" << to_eval_type(tlabel) << "* inteval)";
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
  for(tciter t=taskmgr.first(); t!=taskmgr.last(); ++t) {
    SafePtr<TaskParameters> tparams = t->params();
    const std::string& tlabel = t->label();
    ph_ << define(tlabel,"MAX_STACK_SIZE",tparams->max_stack_size());
    const unsigned int max_vector_stack_size = tparams->max_vector_stack_size();
    if (max_vector_stack_size)
      ph_ << define(tlabel,"MAX_VECTOR_STACK_SIZE",max_vector_stack_size);
    ph_ << define(tlabel,"MAX_HRR_HSRANK",tparams->max_hrr_hsrank());
    ph_ << define(tlabel,"MAX_HRR_LSRANK",tparams->max_hrr_lsrank());
  }

  // libint2_iface.h needs macros to help forming prefixed names in API
  ih_ << "/* Use LIBINT2_PREFIXED_NAME(fncname) to form properly prefixed function name from LIBINT2 API */" << std::endl;
  ih_ << "#define LIBINT2_PREFIXED_NAME(name) __libint2_prefixed_name__(LIBINT2_API_PREFIX,name)" << std::endl;
  ih_ << "#define __libint2_prefixed_name__(prefix,name) __prescanned_prefixed_name__(prefix,name)" << std::endl;
  ih_ << "#define __prescanned_prefixed_name__(prefix,name) prefix##name" << std::endl << std::endl;
  
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
    for(tciter t=taskmgr.first(); t!=taskmgr.last(); ++t,++i) {
      const std::string& tlabel = t->label();

      li_ << lm_decls_[i] << ctext_->open_block();
      li_ << "return " << macro(tlabel,"MAX_STACK_SIZE") << " * VECLEN + "
	  << macro(tlabel,"MAX_VECTOR_STACK_SIZE") << " * VECLEN * ("
	    << macro(tlabel,"MAX_HRR_HSRANK") << " > " << macro(tlabel,"MAX_HRR_LSRANK") << " ? "
	      << macro(tlabel,"MAX_HRR_HSRANK") << " : " << macro(tlabel,"MAX_HRR_LSRANK") << ");" << std::endl;
      li_ << ctext_->close_block();
    }

    i = 0;
    for(tciter t=taskmgr.first(); t!=taskmgr.last(); ++t,++i) {
      const std::string& tlabel = t->label();

      li_ << li_decls_[i] << ctext_->open_block();
      li_ << "if (buf != 0) inteval->stack = buf;" << std::endl << "else ";
      {
	std::string tmp = ctext_->label_to_name(cparams_->api_prefix() + "libint2_need_memory_" + tlabel) + "(max_am)";
	li_ << ctext_->assign("inteval->stack","new LIBINT2_REALTYPE[" + tmp + "]");
      }

      std::string vstack_ptr("inteval->stack + ");
      vstack_ptr += macro(tlabel,"MAX_STACK_SIZE");
      vstack_ptr += " * VECLEN";

      li_ << ctext_->assign("inteval->vstack",vstack_ptr);
      if (cparams_->count_flops()) {
	// set the counter to zero
	li_ << "inteval->nflops = 0;" << endl;
      }
      li_ << ctext_->close_block();
    }

    i = 0;
    for(tciter t=taskmgr.first(); t!=taskmgr.last(); ++t,++i) {
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

