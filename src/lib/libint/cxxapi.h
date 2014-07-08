
#ifndef _libint2_src_lib_libint_simpleapi_h_
#define _libint2_src_lib_libint_simpleapi_h_

#if __cplusplus <= 199711L
# error "Simple Libint API requires C++11 support"
#endif

#include <libint2.h>

namespace libint2 {
  void init() {
    libint2_static_init();
  }
  void cleanup() {
    libint2_static_cleanup();
  }
}

//#include <libint2/engine.h>
#include <engine.h>

#endif /* _libint2_src_lib_libint_simpleapi_h_ */
