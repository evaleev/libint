
#ifndef _libint2_src_bin_libint_smartptr_h_
#define _libint2_src_bin_libint_smartptr_h_

#include <libint2_config.h>

#if HAVE_SHARED_PTR_IN_BOOST
  #include <boost/shared_ptr.hpp>
  #include <boost/enable_shared_from_this.hpp>
  using namespace boost;

  // For now I'll do a cheat since templated typedefs are not standard
  // Should probably at least derive SafePtr from shared_ptr
  #define SafePtr boost::shared_ptr
  #define EnableSafePtrFromThis boost::enable_shared_from_this
  #define SafePtr_from_this shared_from_this
#else
  #include <memory>
  // For now I'll do a cheat since templated typedefs are not standard
  // Should probably at least derive SafePtr from shared_ptr
  #define SafePtr std::shared_ptr
  #define EnableSafePtrFromThis std::enable_shared_from_this
  #define SafePtr_from_this std::shared_from_this
#endif

/** Can be used to determine whether a type is a SafePtr */
template <typename T>
struct IsSafePtr {
  enum { result = false };
};

template <typename T>
struct IsSafePtr< SafePtr<T> > {
  enum { result = true };
};
template <typename T>
struct IsSafePtr< const SafePtr<T> > {
  enum { result = true };
};
template <typename T>
struct IsSafePtr< SafePtr<T>& > {
  enum { result = true };
};
template <typename T>
struct IsSafePtr< const SafePtr<T>& > {
  enum { result = true };
};

#endif

