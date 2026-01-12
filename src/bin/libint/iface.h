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

#ifndef _libint2_src_bin_libint_iface_h_
#define _libint2_src_bin_libint_iface_h_

#include <context.h>
#include <default_params.h>
#include <smart_ptr.h>
#include <task.h>

#include <fstream>
#include <sstream>

namespace libint2 {
/** Libint2Iface is used to generate Libint2 interfaces. The external API
   consists of the main header file libint2.h (not generated), type definition
   file libint2_types.h, parameter file libint2_params.h, function declarations
   in libint2_iface.h, init/cleanup codes libint2_init.cc, libint2_cleanup.cc,
   libint2_static_init.cc, and libint2_static_cleanup.cc. The internal API
   consists of the header file libint2_iface_internal.h, which contains
   prototypes for the top level evaluator functions. The internal API is only
   needed to compile the generated components of Libint2.
  */
class Libint2Iface {
 public:
  typedef std::vector<std::string> Tasks;
  Libint2Iface(const std::shared_ptr<CompilationParameters>& cparams,
               const std::shared_ptr<CodeContext>& ctext);
  ~Libint2Iface();

  /// Writes string s to the types header
  void to_types(const std::string& s) { th_ << s << std::endl; }
  /// Writes string s to the params header
  void to_params(const std::string& s) { ph_ << s << std::endl; }
  /// Writes string s to the iface header
  void to_iface(const std::string& s) { ih_ << s << std::endl; }
  /// Writes string s to the internal iface header
  void to_int_iface(const std::string& s) { ii_ << s << std::endl; }
  /// Writes string s to the static init code
  void to_static_init(const std::string& s) { si_ << s << std::endl; }
  /// Writes string s to the static cleanup code
  void to_static_cleanup(const std::string& s) { sc_ << s << std::endl; }
  /// Writes string s to the Libint_t init code
  // void to_libint_init(const std::string& s) {
  //   li_ << s << endl;
  // }

  std::string macro(const std::string& label) {
    std::string result("LIBINT2_");
    result += label;
    return result;
  }

  std::string macro(const std::string& task_label, const std::string& label) {
    std::string result("LIBINT2_");
    result += label;
    if (task_label != "") {
      result += "_";
      result += task_label;
    }
    return result;
  }

  template <typename T>
  std::string macro_define(const std::string& label, const T& value) {
    oss_.str(null_str_);
    oss_ << "#ifndef " << macro(label) << std::endl;
    oss_ << "# define " << macro(label) << " " << value << std::endl;
    oss_ << "#endif" << std::endl;
    return oss_.str();
  }

  template <typename T>
  std::string macro_define(const std::string& task_label,
                           const std::string& label, const T& value) {
    oss_.str(null_str_);
    oss_ << "#define " << macro(task_label, label) << " " << value << std::endl;
    return oss_.str();
  }

  template <typename T>
  std::string var_declare_v(const std::string& label) {
    return ctext_->declare_v(ctext_->type_name<T>(),
                             ctext_->label_to_name(label), macro("MAX_VECLEN"));
  }

 private:
  std::string null_str_;
  std::ostringstream oss_;
  std::shared_ptr<CompilationParameters> cparams_;
  std::shared_ptr<CodeContext> ctext_;

  // computation-specific functions are libint2_init_xxx, libint2_cleanup_xxx,
  // etc. -- these are their declarations, e.g. "libint2_init_xxx(Libint_t*
  // libint, int max_am, LIBINT2_REALTYPE* buf)"
  std::vector<std::string> li_decls_;  // _init_
  std::vector<std::string> lm_decls_;  // _need_memory_
  std::vector<std::string> lc_decls_;  // _cleanup_

  std::string lf_decl_;  // _init_flopcounter

  typedef std::basic_ofstream<char> fstream;

  fstream th_;
  fstream ph_;
  fstream ih_;
  fstream ii_;
  fstream si_;
  fstream sc_;
  fstream li_;

  void generate_inteval_type(std::ostream& os);
};
};  // namespace libint2

#endif
