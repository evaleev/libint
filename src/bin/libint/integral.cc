
#include <stdexcept>
#include <rr.h>
#include <integral.h>

using namespace std;
using namespace libint2;

const char StaticDefinitions::am_letters[StaticDefinitions::num_am_letters] = "spdfghiklmnoqrtuvwxyz";

namespace libint2 {

  template <>
  bool
  TwoPRep_11_11<CGF>::precomputed() const
  {
    /// (ss|ss)^{(m)} are precomputed 
    if (parent_type::bra_.member(0,0)->zero() && parent_type::bra_.member(1,0)->zero() &&
      parent_type::ket_.member(0,0)->zero() && parent_type::ket_.member(1,0)->zero())
    return true;
    else
      return false;
  }
    
};

