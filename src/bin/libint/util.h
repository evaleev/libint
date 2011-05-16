#ifndef _libint2_src_bin_libint_util_h_
#define _libint2_src_bin_libint_util_h_

#include <string>
#include <stdexcept>
#include <smart_ptr.h>
#include <util_types.h>

namespace libint2 {
  std::string to_string(FunctionPosition pos);
  
  template <class Target, class Source> SafePtr<Target> require_dynamic_cast(const SafePtr<Source>& s) {
    const SafePtr<Target> t = dynamic_pointer_cast<Target,Source>(s);
    if (t == 0)
      throw std::runtime_error("require_dynamic_cast: dynamic case failed");
    return t;
  }
  template <class Target, class Source> const Target* require_dynamic_cast(const Source* s) {
    const Target* t = dynamic_cast<Target*>(s);
    if (t == 0)
      throw std::runtime_error("require_dynamic_cast: dynamic case failed");
    return t;
  }
  
}

#endif /* header guard */
