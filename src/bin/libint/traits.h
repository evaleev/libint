
#include <bfset.h>
#include <smart_ptr.h>
#include <global_macros.h>

#ifndef _libint2_src_bin_libint_traits_h_
#define _libint2_src_bin_libint_traits_h_

namespace libint2 {
  
  template <typename T>
  struct StorageTraits {
    typedef SafePtr<T> StorageType;
    enum { StoredAsPtr = true };
    static const T& const_ref(const StorageType& s) { return *s; };
  };

#if USE_BRAKET_H
  template <>
  struct StorageTraits<CGShell> {
    typedef CGShell StorageType;
    enum { StoredAsPtr = false };
    static const CGShell& const_ref(const StorageType& s) { return s; };
  };
  
  template <>
  struct StorageTraits<CGF> {
    typedef CGF StorageType;
    enum { StoredAsPtr = false };
    static const CGF& const_ref(const StorageType& s) { return s; };
  };
#endif
  
  ///////////
  
  /// Converts Base to a type of the same signature as Ref. For example, if Ref is SafePtr<T> then Base is converted to SafePtr<Base>
  template <typename Ref, typename Base>
  struct ReturnTypeAnalog {
    typedef const Base& result;
  };
  template <typename Ref, typename Base>
  struct ReturnTypeAnalog< SafePtr<Ref>, Base> {
    typedef SafePtr<Base> result;
  };
  
  ///////////
  
  template <typename T>
  struct TypeTraits {
    /// By default, use SafePtr to manage these objects
    typedef typename StorageTraits<T>::StorageType StorageType;
    /// Whether stored as a pointer
    enum { StoredAsPtr = StorageTraits<T>::StoredAsPtr};
    /// Convert an object of StorageType to const T&
    static const T& const_ref(const StorageType& s) { return StorageTraits<T>::const_ref(s); }
  };
};

#endif

