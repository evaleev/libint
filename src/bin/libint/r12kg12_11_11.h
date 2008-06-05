
#ifndef _libint2_src_bin_libint_r12kg121111_h_
#define _libint2_src_bin_libint_r12kg121111_h_

#include <integral.h>
#include <integral_11_11.h>

using namespace std;

namespace libint2 {

  template <>
  inline bool
    GenIntegralSet_11_11<CGF,R12kG12,mType>::this_precomputed() const
    {
#if USE_BRAKET_H
      if (parent_type::bra_.member(0,0).zero() && parent_type::bra_.member(1,0).zero() &&
          parent_type::ket_.member(0,0).zero() && parent_type::ket_.member(1,0).zero())
#else
      if (parent_type::bra_.member(0,0)->zero() && parent_type::bra_.member(1,0)->zero() &&
          parent_type::ket_.member(0,0)->zero() && parent_type::ket_.member(1,0)->zero())
#endif
        return true;
      else
         return false;
    }
  
};

#endif

