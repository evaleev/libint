
#include <flop.h>
#include <rr.h>

using namespace libint2;

namespace libint2 {

  template <>
  unsigned int
  ConvertNumFlops<CGF>(unsigned int nflops) { return nflops; }
  
  template <>
  unsigned int
  ConvertNumFlops<CGShell>(unsigned int nflops) { return 0; }

};
