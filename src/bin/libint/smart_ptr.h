
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#ifndef _libint2_src_bin_libint_smartptr_h_
#define _libint2_src_bin_libint_smartptr_h_

using namespace boost;

// For now I'll do a cheat since templated typedefs are not standard
// Should probably at least derive SafePtr from shared_ptr
#define SafePtr shared_ptr
#define EnableSafePtrFromThis enable_shared_from_this
#define SafePtr_from_this shared_from_this


#endif

