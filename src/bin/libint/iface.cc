
#include <string>
#include <fstream>
#include <iface.h>
#include <default_params.h>

using namespace std;
using namespace libint2;

namespace {
  const char mh_name[] = "libint2.h";
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
                           const SafePtr<CodeContext>& ctext,
                           const Comps& comps) :
  cparams_(cparams), ctext_(ctext), comps_(comps), null_str_(""),
  ph_((cparams_->source_directory() + ph_name).c_str()),
  ih_((cparams_->source_directory() + ih_name).c_str()),
  ii_((cparams_->source_directory() + ii_name).c_str()),
  si_((cparams_->source_directory() + si_name).c_str()),
  sc_((cparams_->source_directory() + sc_name).c_str()),
  li_((cparams_->source_directory() + li_name).c_str())
{
  header_guard_open(ph_,"libint2params");
  header_guard_open(ih_,"libint2iface");
  header_guard_open(ii_,"libint2ifaceint");
  
  ph_ << define("MAX_VECLEN",cparams_->max_vector_length());
  if (cparams_->count_flops())
    ph_ << define("FLOP_COUNT",1);
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

  typedef vector<string> vstype;
  for(vstype::const_iterator ci=comps.begin(); ci <comps.end(); ci++) {
    unsigned int lmax = cparams_->max_am_eri() + 1;
    ostringstream oss;
    oss << "void (*libint2_build_" << *ci
        << "[" << lmax << "][" << lmax << "][" << lmax << "]["
        << lmax << "])(Libint_t *);" << endl;
    ih_ << "extern " << oss.str();
    si_ << oss.str();
  }
  
  // Declare library constructor/destructor
  oss_.str(null_str_);
  oss_ << ctext_->type_name<void>() << " "
       << ctext_->label_to_name("libint2_static_init") << "()";
  std::string si_fdec(oss_.str());
  
  oss_.str(null_str_);
  oss_ << ctext_->type_name<void>() << " "
       << ctext_->label_to_name("libint2_static_cleanup") << "()";
  std::string sc_fdec(oss_.str());

  ih_ << si_fdec << ctext_->end_of_stat() << endl;
  ih_ << sc_fdec << ctext_->end_of_stat() << endl;

  // Declare Libint_t constructor/destructor (specific to the type of computation)
  for(vstype::const_iterator ci=comps.begin(); ci <comps.end(); ci++) {
    oss_.str(null_str_);
    oss_ << ctext_->type_name<void>() << " "
	 << ctext_->label_to_name("libint2_init_" + *ci) << "(Libint_t* libint, int max_am, LIBINT2_REALTYPE* buf)";
    std::string li_fdec(oss_.str());
    li_decls_.push_back(li_fdec);
  
    oss_.str(null_str_);
    oss_ << ctext_->type_name<size_t>() << " "
	 << ctext_->label_to_name("libint2_need_memory_" + *ci) << "(int max_am)";
    std::string lm_fdec(oss_.str());
    lm_decls_.push_back(lm_fdec);

    oss_.str(null_str_);
    oss_ << ctext_->type_name<void>() << " "
	 << ctext_->label_to_name("libint2_cleanup_" + *ci) << "(Libint_t* libint)";
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
  LibraryParameters& lparams = LibraryParameters::get_library_params();
  ph_ << define("MAX_STACK_SIZE",lparams.max_stack_size());
  const unsigned int max_vector_stack_size = lparams.max_vector_stack_size();
  if (max_vector_stack_size)
    ph_ << define("MAX_VECTOR_STACK_SIZE",max_vector_stack_size);
  ph_ << define("MAX_HRR_HSRANK",lparams.max_hrr_hsrank());
  ph_ << define("MAX_HRR_LSRANK",lparams.max_hrr_lsrank());
  
  header_guard_close(ph_);
  header_guard_close(ih_);
  header_guard_close(ii_);

  std::ostringstream oss;
  oss << ctext_->close_block() << ctext_->code_postfix() << endl << endl;
  std::string pfix = oss.str();

  si_ << pfix;
  sc_ << pfix;

  // Define Libint_t constructor/destructor (specific to the type of computation)
  const unsigned int ncomps = comps_.size();
  for(unsigned int i=0; i<ncomps; ++i) {
    li_ << lm_decls_[i] << ctext_->open_block();
    li_ << "return LIBINT2_MAX_STACK_SIZE * VECLEN + LIBINT2_MAX_VECTOR_STACK_SIZE * VECLEN * (LIBINT2_MAX_HRR_HSRANK > LIBINT2_MAX_HRR_LSRANK ? LIBINT2_MAX_HRR_HSRANK : LIBINT2_MAX_HRR_LSRANK);" << std::endl;
    li_ << ctext_->close_block();
  }
  for(unsigned int i=0; i<ncomps; ++i) {
    li_ << li_decls_[i] << ctext_->open_block();
    li_ << "if (buf != 0) libint->stack = buf;" << std::endl << "else ";
    {
      std::string tmp = ctext_->label_to_name("libint2_need_memory_" + comps_[i]) + "(max_am)";
      li_ << ctext_->assign("libint->stack","new LIBINT2_REALTYPE[" + tmp + "]");
    }
    li_ << ctext_->assign("libint->vstack","libint->stack + LIBINT2_MAX_STACK_SIZE * VECLEN");
    if (cparams_->count_flops()) {
      // set the counter to zero
      li_ << "libint->nflops = 0;" << endl;
    }
    li_ << ctext_->close_block();
  }
  for(unsigned int i=0; i<ncomps; ++i) {
    li_ << lc_decls_[i] << ctext_->open_block();
    li_ << "delete[] libint->stack;" << std::endl;
    li_ << ctext_->assign("libint->stack","0");
    li_ << ctext_->assign("libint->vstack","0");
    if (cparams_->count_flops()) {
      // set the counter to zero
      li_ << "libint->nflops = 0;" << endl;
    }
    li_ << ctext_->close_block();
  }

  li_ << ctext_->code_postfix() << std::endl << std::endl;

  ph_.close();
  ih_.close();
  ii_.close();
  si_.close();
  sc_.close();
  li_.close();
}

