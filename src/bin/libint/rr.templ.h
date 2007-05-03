
#ifndef _libint2_src_bin_libint_rrtemplh_h_
#define _libint2_src_bin_libint_rrtemplh_h_

#include <rr.h>

using namespace libint2;

template <class RR> bool
RecurrenceRelation::register_with_rrstack() const {
  // only register RRs with for shell sets
  if (TrivialBFSet<typename RR::BasisFunctionType>::result)
    return false;
  SafePtr<RRStack> rrstack = RRStack::Instance();
  SafePtr<RR> this_ptr =
    const_pointer_cast<RR,const RR>(
      static_pointer_cast<const RR, const RecurrenceRelation>(
        EnableSafePtrFromThis<RecurrenceRelation>::SafePtr_from_this()
      )
    );
  rrstack->find(this_ptr);
  return true;
}

#endif
