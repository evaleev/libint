
#include <iostream>
#include <stdexcept>
#include <integral.h>

using namespace std;
using namespace libint2;

namespace libint2 {

  template <>
  bool
  TwoPRep_11_11<CGShell>::this_precomputed() const
  {
    return false;
  }
    
  template <>
  bool
  TwoPRep_11_11<CGF>::this_precomputed() const
  {
    /// (ss|ss)^{(m)} are precomputed 
    if (parent_type::bra_.member(0,0)->zero() && parent_type::bra_.member(1,0)->zero() &&
      parent_type::ket_.member(0,0)->zero() && parent_type::ket_.member(1,0)->zero())
      return true;
    else
      return false;
  }

};

