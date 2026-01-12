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

#include <default_params.h>
#include <iface.h>

#include <fstream>
#include <string>

using namespace std;
using namespace libint2;

namespace {
const char th_name[] = "libint2_types.h";
const char ph_name[] = "libint2_params.h";
const char ih_name[] = "libint2_iface.h";
const char ii_name[] = "libint2_iface_internal.h";
const char si_name[] = "libint2_static_init.cc";
const char sc_name[] = "libint2_static_cleanup.cc";
const char li_name[] = "libint2_init.cc";

inline void header_guard_open(std::ostream& os, const std::string& label) {
  os << "#ifndef _libint2_" << label << "_h_" << endl
     << "#define _libint2_" << label << "_h_" << endl
     << endl;
}
inline void header_guard_close(std::ostream& os) {
  os << "#endif" << endl << endl;
}

};  // namespace

Libint2Iface::Libint2Iface(
    const std::shared_ptr<CompilationParameters>& cparams,
    const std::shared_ptr<CodeContext>& ctext)
    : null_str_(""),
      oss_(),
      cparams_(cparams),
      ctext_(ctext),
      th_((cparams_->source_directory() + th_name).c_str()),
      ph_((cparams_->source_directory() + ph_name).c_str()),
      ih_((cparams_->source_directory() + ih_name).c_str()),
      ii_((cparams_->source_directory() + ii_name).c_str()),
      si_((cparams_->source_directory() + si_name).c_str()),
      sc_((cparams_->source_directory() + sc_name).c_str()),
      li_((cparams_->source_directory() + li_name).c_str()) {
  th_ << ctext_->copyright();
  ph_ << ctext_->copyright();
  ih_ << ctext_->copyright();
  ii_ << ctext_->copyright();
  si_ << ctext_->copyright();
  sc_ << ctext_->copyright();
  li_ << ctext_->copyright();

  header_guard_open(th_, "libint2types");
  header_guard_open(ph_, "libint2params");
  header_guard_open(ih_, "libint2iface");
  header_guard_open(ii_, "libint2ifaceint");

  ph_ << macro_define("API_PREFIX", cparams_->api_prefix());
  ph_ << macro_define("MAX_VECLEN", cparams_->max_vector_length());
  ph_ << macro_define("ALIGN_SIZE", cparams_->align_size());
  if (cparams_->count_flops()) ph_ << macro_define("FLOP_COUNT", 1);
  if (cparams_->profile()) ph_ << macro_define("PROFILE", 1);
  if (cparams_->accumulate_targets()) ph_ << macro_define("ACCUM_INTS", 1);
  const std::string realtype(cparams_->realtype());
  {
    // does LIBINT_USER_DEFINED_REAL need extra include statements?
#ifdef LIBINT_USER_DEFINED_REAL_INCLUDES
    ph_ << LIBINT_USER_DEFINED_REAL_INCLUDES << std::endl;
#endif
    ph_ << macro_define("REALTYPE", realtype);
  }
  if (cparams_->contracted_targets()) ph_ << macro_define("CONTRACTED_INTS", 1);

  ih_ << "#ifdef __cplusplus\n# include <cstddef>\n#else\n# include "
         "<stddef.h>\n#endif"
      << endl
      << ctext_->code_prefix();

  oss_.str(null_str_);
  oss_ << ctext_->std_header() << "#include <" << ih_name << ">" << endl
       << "#include <" << ii_name << ">" << endl
       << "#include <cstddef>" << endl
       << "#include <cassert>" << endl
       << "#include <cstdlib>" << endl
       << ctext_->code_prefix();
  std::string pfix = oss_.str();
  si_ << pfix;
  sc_ << pfix;
  li_ << "#include <libint2/util/memory.h>" << endl
      << pfix;  // moved from libint2.h to here

  // print out declarations for the array of pointers to evaluator functions
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  typedef LibraryTaskManager::TasksCIter tciter;
  for (tciter t = taskmgr.first(); t != taskmgr.plast(); ++t) {
    const std::string& tlabel = t->label();
    const unsigned int nbf = cparams_->num_bf(tlabel);

    ostringstream oss;
    oss << "void (*" << ctext->label_to_name(cparams->api_prefix())
        << "libint2_build_" << tlabel;
    for (unsigned int c = 0; c < nbf; ++c) {
      const unsigned int lmax =
          const_cast<const CompilationParameters*>(cparams_.get())
              ->max_am(tlabel, c);
      oss << "[" << lmax + 1 << "]";
      std::cout << "task=" << tlabel << " center=" << c << " lmax=" << lmax
                << std::endl;
    }
    oss << ")(" << ctext_->const_modifier() << ctext_->inteval_type_name(tlabel)
        << "*);" << endl;
    ih_ << "extern " << oss.str();
    si_ << oss.str();
  }

  // Declare library constructor/destructor
  oss_.str(null_str_);
  oss_ << ctext_->type_name<void>() << " "
       << ctext_->label_to_name(cparams->api_prefix() + "libint2_static_init")
       << "()";
  std::string si_fdec(oss_.str());

  oss_.str(null_str_);
  oss_ << ctext_->type_name<void>() << " "
       << ctext_->label_to_name(cparams->api_prefix() +
                                "libint2_static_cleanup")
       << "()";
  std::string sc_fdec(oss_.str());

  ih_ << si_fdec << ctext_->end_of_stat() << endl;
  ih_ << sc_fdec << ctext_->end_of_stat() << endl;

  // Declare evaluator constructor/destructor (specific to the type of
  // computation)
  for (tciter t = taskmgr.first(); t != taskmgr.plast(); ++t) {
    const std::string& tlabel = t->label();

    oss_.str(null_str_);
    oss_ << ctext_->type_name<void>() << " "
         << ctext_->label_to_name(cparams->api_prefix() + "libint2_init_" +
                                  tlabel)
         << "(" << ctext_->inteval_type_name(tlabel)
         << "* inteval, int max_am, void* buf)";
    std::string li_fdec(oss_.str());
    li_decls_.push_back(li_fdec);

    oss_.str(null_str_);
    oss_ << ctext_->type_name<size_t>() << " "
         << ctext_->label_to_name(cparams->api_prefix() +
                                  "libint2_need_memory_" + tlabel)
         << "(int max_am)";
    std::string lm_fdec(oss_.str());
    lm_decls_.push_back(lm_fdec);

    oss_.str(null_str_);
    oss_ << ctext_->type_name<void>() << " "
         << ctext_->label_to_name(cparams->api_prefix() + "libint2_cleanup_" +
                                  tlabel)
         << "(" << ctext_->inteval_type_name(tlabel) << "* inteval)";
    std::string lc_fdec(oss_.str());
    lc_decls_.push_back(lc_fdec);

    ih_ << li_fdec << ctext_->end_of_stat() << endl;
    ih_ << lm_fdec << ctext_->end_of_stat() << endl;
    ih_ << lc_fdec << ctext_->end_of_stat() << endl;
  }

  // if counting flops, need additional initializer function
  if (cparams_->count_flops()) {
    oss_.str(null_str_);
    oss_ << "#ifdef __cplusplus\n#ifdef LIBINT2_FLOP_COUNT\nextern \"C++\" "
            "template <typename EvalType> void "
         << ctext_->label_to_name(cparams->api_prefix() +
                                  "libint2_init_flopcounter")
         << "(EvalType* inteval_vector, int inteval_vector_size)"
         << ctext_->open_block();
    // TODO convert to ForLoop object
    oss_ << "for(int v=1; v!=inteval_vector_size; ++v)" << ctext_->open_block()
         << ctext_->assign("inteval_vector[v].nflops",
                           "inteval_vector[0].nflops")
         << ctext_->close_block() << ctext_->close_block();
    oss_ << "#endif\n#endif\n";

    lf_decl_ = oss_.str();
    ih_ << lf_decl_ << ctext_->end_of_stat() << endl;
  }

  ih_ << ctext_->code_postfix() << endl;

  si_ << si_fdec << ctext_->open_block();
  sc_ << sc_fdec << ctext_->open_block();
}

Libint2Iface::~Libint2Iface() {
  // For each task, print out defines for stack dimensions
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  typedef LibraryTaskManager::TasksCIter tciter;
  for (tciter t = taskmgr.first(); t != taskmgr.plast(); ++t) {
    std::shared_ptr<TaskParameters> tparams = t->params();
    const std::string& tlabel = t->label();
    ph_ << macro_define(tlabel, "NUM_TARGETS", tparams->max_ntarget());
    const unsigned int max_am = tparams->max_am();
    for (unsigned int am = 0; am <= max_am; ++am) {
      {
        std::ostringstream oss;
        oss << "MAX_STACK_SIZE_" << am;
        ph_ << macro_define(tlabel, oss.str(), tparams->max_stack_size(am));
      }
      {
        std::ostringstream oss;
        oss << "MAX_VECTOR_STACK_SIZE_" << am;
        const unsigned int max_vector_stack_size =
            tparams->max_vector_stack_size(am);
        ph_ << macro_define(tlabel, oss.str(), max_vector_stack_size);
      }
      {
        std::ostringstream oss;
        oss << "MAX_HRR_HSRANK_" << am;
        ph_ << macro_define(tlabel, oss.str(), tparams->max_hrr_hsrank(am));
      }
      {
        std::ostringstream oss;
        oss << "MAX_HRR_LSRANK_" << am;
        ph_ << macro_define(tlabel, oss.str(), tparams->max_hrr_lsrank(am));
      }
    }
  }

  // For each task, generate the evaluator type
  th_ << "#include <libint2/util/vector.h>" << std::endl;
  th_ << "#include <libint2/util/intrinsic_types.h>" << std::endl;
  th_ << "#include <libint2/util/intrinsic_operations.h>" << std::endl;
  th_ << "#include <libint2/util/timer.h>"
      << std::endl;  // in case LIBINT2_PROFILE is on
  generate_inteval_type(th_);

  // libint2_iface.h needs macros to help forming prefixed names in API
  ih_ << ctext_->comment(
             "Use LIBINT2_PREFIXED_NAME(fncname) to form properly prefixed "
             "function name from LIBINT2 API")
      << std::endl;
  ih_ << "#define LIBINT2_PREFIXED_NAME(name) "
         "__libint2_prefixed_name__(LIBINT2_API_PREFIX,name)"
      << std::endl;
  ih_ << "#define __libint2_prefixed_name__(prefix,name) "
         "__prescanned_prefixed_name__(prefix,name)"
      << std::endl;
  ih_ << "#define __prescanned_prefixed_name__(prefix,name) prefix##name"
      << std::endl;
  // also need macros to test what's in the evaluator
  ih_ << ctext_->comment(
             "Use LIBINT2_PREFIXED_NAME(fncname) to form properly prefixed "
             "function name from LIBINT2 API")
      << std::endl;
  ih_ << "#define LIBINT2_DEFINED(taskname,symbol) "
         "__prescanned_libint2_defined__(taskname,symbol)"
      << std::endl;
  if (cparams_->single_evaltype())
    ih_ << "#define __prescanned_libint2_defined__(taskname,symbol) "
           "LIBINT2_DEFINED_##symbol"
        << std::endl
        << std::endl;
  else
    ih_ << "#define __prescanned_libint2_defined__(taskname,symbol) "
           "LIBINT2_DEFINED_##symbol##_##taskname"
        << std::endl
        << std::endl;

  header_guard_close(th_);
  header_guard_close(ph_);
  header_guard_close(ih_);
  header_guard_close(ii_);

  std::ostringstream oss;
  oss << ctext_->close_block() << ctext_->code_postfix() << endl << endl;
  std::string pfix = oss.str();

  si_ << pfix;
  sc_ << pfix;

  // Define Libint_t constructor/destructor (specific to the type of
  // computation)
  {
    unsigned int i = 0;
    for (tciter t = taskmgr.first(); t != taskmgr.plast(); ++t, ++i) {
      const std::string& tlabel = t->label();

      li_ << lm_decls_[i] << ctext_->open_block();
      const unsigned int max_am = t->params()->max_am();
      for (unsigned int am = 0; am <= max_am; ++am) {
        std::string ss, vss, hsr, lsr;
        {
          std::ostringstream oss;
          oss << "MAX_STACK_SIZE_" << am;
          ss = oss.str();
        }
        {
          std::ostringstream oss;
          oss << "MAX_VECTOR_STACK_SIZE_" << am;
          vss = oss.str();
        }
        {
          std::ostringstream oss;
          oss << "MAX_HRR_HSRANK_" << am;
          hsr = oss.str();
        }
        {
          std::ostringstream oss;
          oss << "MAX_HRR_LSRANK_" << am;
          lsr = oss.str();
        }

        li_ << "assert(max_am <= " << max_am << ");" << std::endl;

        li_ << "if (max_am == " << am << ") return " << macro(tlabel, ss)
            << " * " << macro("MAX_VECLEN") << " + " << macro(tlabel, vss)
            << " * " << macro("MAX_VECLEN") << " * (" << macro(tlabel, hsr)
            << " > " << macro(tlabel, lsr) << " ? " << macro(tlabel, hsr)
            << " : " << macro(tlabel, lsr) << ");" << std::endl;
      }
      li_ << "return 0; // unreachable" << std::endl;
      li_ << ctext_->close_block();
    }

    i = 0;
    for (tciter t = taskmgr.first(); t != taskmgr.plast(); ++t, ++i) {
      const std::string& tlabel = t->label();

      li_ << li_decls_[i] << ctext_->open_block();
      li_ << "if (buf != 0) inteval->stack = "
             "reinterpret_cast<LIBINT2_REALTYPE*>(buf);"
          << std::endl
          << "else " << std::endl;
      {
        std::string tmp =
            ctext_->label_to_name(cparams_->api_prefix() +
                                  "libint2_need_memory_" + tlabel) +
            "(max_am)";

        // no posix_memalign? use new, with default alignment
        li_ << ctext_->assign(
            "inteval->stack",
            std::string("libint2::malloc<LIBINT2_REALTYPE>(") + tmp +
                std::string(")"));
      }

      const unsigned int max_am = t->params()->max_am();
      for (unsigned int am = 0; am <= max_am; ++am) {
        std::string ss;
        {
          std::ostringstream oss;
          oss << "MAX_STACK_SIZE_" << am;
          ss = oss.str();
        }

        li_ << "assert(max_am <= " << max_am << ");" << std::endl;
        li_ << "if (max_am == " << am << ")" << std::endl;
        std::string vstack_ptr("inteval->stack + ");
        vstack_ptr += macro(tlabel, ss);
        vstack_ptr += " * ";
        vstack_ptr += macro("MAX_VECLEN");
        li_ << ctext_->assign("inteval->vstack", vstack_ptr);
      }

      if (cparams_->count_flops()) {
        // allocate the counter and set it to zero
        li_ << "inteval->nflops = new " << macro("UINT_LEAST64") << ";" << endl;
        li_ << "inteval->nflops[0] = 0;" << endl;
      }
      if (cparams_->profile()) {  // zero out the timers
        li_ << ctext_->macro_if("LIBINT2_CPLUSPLUS_STD >= 2011");
        li_ << "inteval->timers = new libint2::Timers<2>;" << endl;
        li_ << "inteval->timers->clear();" << endl;
        li_ << ctext_->macro_endif();  // >= C++11
      }
      li_ << ctext_->close_block();
    }

    i = 0;
    for (tciter t = taskmgr.first(); t != taskmgr.plast(); ++t, ++i) {
      li_ << lc_decls_[i] << ctext_->open_block();

      li_ << "free(inteval->stack);\n";
      li_ << ctext_->assign("inteval->stack", "0");
      li_ << ctext_->assign("inteval->vstack", "0");
      if (cparams_->count_flops()) {
        // free the counter and set the pointer to zero
        li_ << "delete inteval->nflops;" << endl;
        li_ << "inteval->nflops = 0;" << endl;
      }
      if (cparams_->profile()) {
        li_ << ctext_->macro_if("LIBINT2_CPLUSPLUS_STD >= 2011");
        li_ << "delete inteval->timers;" << endl;
        li_ << "inteval->timers = 0;" << endl;
        li_ << ctext_->macro_endif();  // >= C++11
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

namespace libint2 {

/// Parses the symbol if it is composed of a prefix followed by a number.
class Parser_prefixN {
 public:
  Parser_prefixN(const std::string& symbol) : match_(false), N_(0) {
    std::string N_str;
    bool minus = false;

    auto prefix_rbegin = symbol.rbegin();
    for (auto ri = symbol.rbegin(); ri != symbol.rend(); ++ri) {
      if (isdigit(*ri)) {
        match_ = true;
        N_str.insert(0, 1, *ri);
      } else {  // no more digits? check for a minus sign
        prefix_rbegin = ri;
        if (*ri == '_') {
          std::string token("_");
          ++ri;
          for (int i = 0; i < 6 && ri != symbol.rend(); ++ri, ++i) {
            token.push_back(*ri);
          }
          if (token == "_sunim_") {  // minus reversed
            minus = true;
            prefix_rbegin = ri;
          }
        }
        break;
      }
    }

    if (match_) {
      N_ = atoi(N_str.c_str());
      if (minus) N_ *= -1;

      prefix_ = std::string(prefix_rbegin, symbol.rend());
      std::reverse(prefix_.begin(), prefix_.end());
    }
  }

  /// returns true if the pattern matched
  bool match() const { return match_; }
  /// returns N (only valid if match()==true)
  int N() const { return N_; }
  /// returns prefix (only valid if match()==true)
  const std::string& prefix() const { return prefix_; }

 private:
  bool match_;
  std::string prefix_;
  int N_;
};
}  // namespace libint2

void Libint2Iface::generate_inteval_type(std::ostream& os) {
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();

  //
  // If need to generate single type for all tasks, take a union of all symbols
  // else process each task separately
  //
  typedef LibraryTaskManager::TasksCIter tciter;
  const tciter tend =
      cparams_->single_evaltype() ? taskmgr.first() + 1 : taskmgr.plast();
  for (tciter t = taskmgr.first(); t != tend; ++t) {
    const std::shared_ptr<TaskExternSymbols> tsymbols = t->symbols();

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
      for (tciter t = taskmgr.first(); t != tend; ++t) {
        const std::shared_ptr<TaskExternSymbols> tsymbols = t->symbols();
        composite_symbols.add(tsymbols->symbols());
      }
      symbols = composite_symbols.symbols();
      tlabel = "";
    } else {
      symbols = tsymbols->symbols();
      tlabel = t->label();
    }

    // Spend some effort to ensure reasonable ordering/grouping of symbols in
    // the list, e.g. all symbols of form PrefixIndex (e.g. crazyprefixN) will
    // be grouped together in the order of increasing Index.
    SymbolList ordered_symbols;
    // 1) convert all symbols to valid code
    // 2) then scan for symbols matching pattern prefixN, for each prefix
    // determine range of N
    std::map<std::string, std::pair<int, int> > prefix_symbols;
    for (auto s : symbols) {
      Parser_prefixN parser(ctext_->label_to_name(s));
      if (parser.match()) {
        if (prefix_symbols.find(parser.prefix()) == prefix_symbols.end()) {
          prefix_symbols[parser.prefix()] =
              std::make_pair(parser.N(), parser.N());
        } else {
          const int N = parser.N();
          auto Nlimits = prefix_symbols[parser.prefix()];
          if (N < Nlimits.first) Nlimits.first = N;
          if (N > Nlimits.second) Nlimits.second = N;
          prefix_symbols[parser.prefix()] = Nlimits;
        }
      }
    }

    // 3) for each prefix iterate over the corresponding range, and add the
    // matching symbols in the same order
    for (auto pi = prefix_symbols.begin(); pi != prefix_symbols.end(); ++pi) {
      auto prefix = pi->first;
      auto Nlimits = pi->second;
      for (int N = Nlimits.first; N <= Nlimits.second; ++N) {
        std::ostringstream oss;
        oss << prefix << (N < 0 ? "-" : "") << abs(N);
        std::string code_symbol = ctext_->label_to_name(oss.str());
        for (auto v : symbols) {
          if (code_symbol == ctext_->label_to_name(v)) {
            ordered_symbols.push_back(code_symbol);
          }
        }
      }
    }
    // 4) the rest of symbols can appear in any order
    //    in practice, the order is lexicographic, hence xyz components are
    //    automatically ordered x,y,z
    for (auto s : symbols) {
      if (std::find(ordered_symbols.begin(), ordered_symbols.end(),
                    ctext_->label_to_name(s)) == ordered_symbols.end())
        ordered_symbols.push_back(s);
    }

    // now dump all symbols into the evaluator data type definition
    for (auto s = ordered_symbols.begin(); s != ordered_symbols.end(); ++s) {
      // for each extrnal symbol #define a macro
      std::string tmplabel("DEFINED_");
      tmplabel += ctext_->label_to_name(*s);
      os << macro_define(tlabel, tmplabel, 1);
      // all symbols are doubles
      os << var_declare_v<double>(*s);
    }

    // Declare members common to all evaluators
    os << ctext_->comment("Scratch buffer to hold intermediates") << std::endl;
    os << ctext_->declare(
        ctext_->mutable_modifier() + ctext_->type_name<double*>(),
        std::string("stack"));

    os << ctext_->comment(
              "Buffer to hold vector intermediates. Only used by set-level RR "
              "code if it is vectorized linewise")
       << std::endl;
    os << ctext_->declare(
        ctext_->mutable_modifier() + ctext_->type_name<double*>(),
        std::string("vstack"));

    os << ctext_->comment(
              "On completion, this contains pointers to computed targets")
       << std::endl;
    if (cparams_->single_evaltype()) {
      // figure out the maximum number of targets
      unsigned int max_ntargets = 0;
      const tciter tend = taskmgr.plast();
      for (tciter t = taskmgr.first(); t != tend; ++t) {
        std::shared_ptr<TaskParameters> tparams = t->params();
        max_ntargets = std::max(max_ntargets, tparams->max_ntarget());
      }
      ostringstream oss;
      oss << max_ntargets;
      os << ctext_->declare_v(
          ctext_->mutable_modifier() + ctext_->type_name<double*>(),
          std::string("targets"), oss.str());
    } else {
      os << ctext_->declare_v(
          ctext_->mutable_modifier() + ctext_->type_name<double*>(),
          std::string("targets"), macro(tlabel, "NUM_TARGETS"));
    }

    os << ctext_->comment(
              "Actual vector length. Not to exceed MAX_VECLEN! If MAX_VECLEN "
              "is 1 then veclen is not used")
       << std::endl;
    os << ctext_->declare(ctext_->type_name<int>(), std::string("veclen"));

    os << ctext_->macro_if(macro("FLOP_COUNT"));
    os << ctext_->comment(
              "FLOP counter. Libint must be configured with "
              "--enable-flop-counter to allow FLOP counting. It is user's "
              "reponsibility to set zero nflops before computing integrals.")
       << std::endl;
    os << ctext_->declare(ctext_->mutable_modifier() + macro("UINT_LEAST64*"),
                          std::string("nflops"));
    os << ctext_->macro_endif();

    os << ctext_->macro_if(macro("PROFILE"));
    os << ctext_->macro_if("LIBINT2_CPLUSPLUS_STD >= 2011");
    os << ctext_->comment(
              "profiling timers. Libint must be configured with "
              "--enable-profile to allow profiling.")
       << std::endl;
    os << "#ifdef __cplusplus" << std::endl
       << ctext_->declare(
              ctext_->mutable_modifier() + "libint2::Timers<2>*",
              std::string("timers"))  // 1 timer for HRR and 1 timer for VRR
       << "#else // timers are not accessible from C" << std::endl
       << "  void* timers;" << std::endl
       << "#endif" << std::endl;
    os << ctext_->macro_endif();  // >= C++11
    os << ctext_->macro_endif();

    os << ctext_->macro_if(macro("ACCUM_INTS"));
    os << ctext_->comment(
              "If libint was configured with --enable-accum-ints then the "
              "target integrals are accumulated. To zero out the targets "
              "automatically before the computation, set this to nonzero.")
       << std::endl;
    os << ctext_->declare(ctext_->type_name<int>(),
                          std::string("zero_out_targets"));
    os << ctext_->macro_endif();

    os << ctext_->macro_if(macro("CONTRACTED_INTS"));
    os << ctext_->comment(
              "If libint was configured with --enable-contracted-ints then "
              "contracted integrals are supported. Set this parameter to the "
              "total number of primitive combinations.")
       << std::endl;
    os << ctext_->declare(ctext_->type_name<int>(), std::string("contrdepth"));
    os << ctext_->macro_endif();

    // Epilogue
    os << "} " << ctext_->inteval_type_name(tlabel) << ctext_->end_of_stat()
       << std::endl;
  }

  // If generating single evaluator type, create aliases from the specialized
  // types to the actual type
  if (cparams_->single_evaltype()) {
    const tciter tend = taskmgr.plast();
    for (tciter t = taskmgr.first(); t != tend; ++t) {
      const std::string& tlabel = t->label();
      os << "typedef " << ctext_->inteval_gen_type_name() << " "
         << ctext_->inteval_spec_type_name(tlabel) << ctext_->end_of_stat()
         << std::endl;
    }
  }
}
