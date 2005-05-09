
#include <ostream>
#include <sstream>
#include <smart_ptr.h>
#include <default_params.h>
#include <context.h>

#ifndef _libint2_src_bin_libint_iface_h_
#define _libint2_src_bin_libint_iface_h_

using namespace std;

namespace libint2 {
  /** Libint2Iface is used to generate Libint2 interfaces. The external API consists of 
      the header file libint2.h (it is not generated yet), parameter file
      libint2_params.h, function declarations in libint2_iface.h,
      init/cleanup codes libint2_init.cc, libint2_cleanup.cc, libint2_static_init.cc,
      and libint2_static_cleanup.cc. The internal API consists of the header
      file libint2_iface_internal.h, which contains prototypes for the top level
      evaluator functions. The internal API is only needed to compile the generated
      components of Libint2.
    */
    class Libint2Iface {
      public:
      typedef std::vector<std::string> Comps;
      Libint2Iface(const SafePtr<CompilationParameters>& cparams,
                   const SafePtr<CodeContext>& ctext,
                   const Comps& comps);
      ~Libint2Iface();
      
      /// Writes string s to the params header
      void to_params(const std::string& s) {
        ph_ << s << endl;
      }
      /// Writes string s to the iface header
      void to_iface(const std::string& s) {
        ih_ << s << endl;
      }
      /// Writes string s to the internal iface header
      void to_int_iface(const std::string& s) {
        ii_ << s << endl;
      }
      /// Writes string s to the static init code
      void to_static_init(const std::string& s) {
        si_ << s << endl;
      }
      /// Writes string s to the static cleanup code
      void to_static_cleanup(const std::string& s) {
        sc_ << s << endl;
      }
      /// Writes string s to the Libint_t init code
      void to_libint_init(const std::string& s) {
        li_ << s << endl;
      }
      /// Writes string s to the Libint_t cleanup code
      void to_libint_cleanup(const std::string& s) {
        lc_ << s << endl;
      }
      
      template <typename T> const std::string define(const std::string& label, const T& value) {
        oss_ .str(null_str_);
        oss_ << "#define LIBINT2_" << label << " " << value << endl;
        return oss_.str();
      }
      
      private:
      std::string null_str_;
      std::ostringstream oss_;
      SafePtr<CompilationParameters> cparams_;
      SafePtr<CodeContext> ctext_;
      Comps comps_;
      
      typedef std::basic_ofstream<char> fstream;
      
      fstream ph_;
      fstream ih_;
      fstream ii_;
      fstream si_;
      fstream sc_;
      fstream li_;
      fstream lc_;
      
    };
};

#endif

