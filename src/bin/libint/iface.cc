
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
  const char lc_name[] = "libint2_cleanup.cc";
  
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
  li_((cparams_->source_directory() + li_name).c_str()),
  lc_((cparams_->source_directory() + lc_name).c_str())
{
  header_guard_open(ph_,"libint2params");
  header_guard_open(ih_,"libint2iface");
  header_guard_open(ii_,"libint2ifaceint");
  
  ph_ << define("MAX_VECLEN",cparams_->max_vector_length());
  
  ih_ << ctext_->code_prefix();

  oss_.str(null_str_);
  oss_ << ctext_->std_header() << "#include <" << ih_name << ">" << endl
                               << "#include <" << ii_name << ">" << endl
                               << ctext_->code_prefix();
  std::string pfix = oss_.str();
  si_ << pfix;
  sc_ << pfix;
  li_ << pfix;
  lc_ << pfix;

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
  
  oss_.str(null_str_);
  oss_ << ctext_->type_name<void>() << " "
       << ctext_->label_to_name("libint2_static_init") << "()";
  std::string si_fdec(oss_.str());
  
  oss_.str(null_str_);
  oss_ << ctext_->type_name<void>() << " "
       << ctext_->label_to_name("libint2_static_cleanup") << "()";
  std::string sc_fdec(oss_.str());
  
  oss_.str(null_str_);
  oss_ << ctext_->type_name<void>() << " "
       << ctext_->label_to_name("libint2_init") << "(Libint_t* libint, int max_am)";
  std::string li_fdec(oss_.str());
  
  oss_.str(null_str_);
  oss_ << ctext_->type_name<void>() << " "
       << ctext_->label_to_name("libint2_cleanup") << "(Libint_t* libint)";
  std::string lc_fdec(oss_.str());
  
  ih_ << si_fdec << ctext_->end_of_stat() << endl;
  ih_ << sc_fdec << ctext_->end_of_stat() << endl;
  ih_ << li_fdec << ctext_->end_of_stat() << endl;
  ih_ << lc_fdec << ctext_->end_of_stat() << endl;
  ih_ << ctext_->code_postfix() << endl;
  
  si_ << si_fdec << ctext_->open_block();
  sc_ << sc_fdec << ctext_->open_block();
  li_ << li_fdec << ctext_->open_block();
  lc_ << lc_fdec << ctext_->open_block();
  
#if UPDATE_FLOP_COUNTER
  // set the counter to zero
  li_ << "libint->nflops = 0;" << endl;
#endif
}

Libint2Iface::~Libint2Iface()
{
  LibraryParameters& lparams = LibraryParameters::get_library_params();
  ph_ << define("MAX_STACK_SIZE",lparams.max_stack_size());
  
  header_guard_close(ph_);
  header_guard_close(ih_);
  header_guard_close(ii_);
  
  std::ostringstream oss;
  oss << ctext_->close_block() << ctext_->code_postfix() << endl << endl;
  std::string pfix = oss.str();
  si_ << pfix;
  sc_ << pfix;
  li_ << pfix;
  lc_ << pfix;
  
  ph_.close();
  ih_.close();
  ii_.close();
  si_.close();
  sc_.close();
  li_.close();
  lc_.close();
}

