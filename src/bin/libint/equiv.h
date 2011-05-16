
#include <smart_ptr.h>
#include <global_macros.h>

#ifndef _libint2_src_bin_libint_equiv_h_
#define _libint2_src_bin_libint_equiv_h_

namespace libint2 {

  /**
  PtrEquiv<T> provides a set of comparison functions named 'equiv' which take
  as arguments a mix of references, regular pointers, and smart pointers to T
  and it's various expected relatives. T must define the type of its parent
  publicly as 'parent_type'.
  */
  template <class T>
  class PtrEquiv {
  
  public:

    /// A shortcut for T::parent_type
    typedef typename T::parent_type P;
    
    static bool equiv(const T& a, const T& b) {
      return a==b;
    }
    
    static bool equiv(const SafePtr<T>& a, const SafePtr<T>& b) {
      return a->operator==(*b.get());
    }

    static bool equiv(const T* a, const SafePtr<T>& b) {
      return a->operator==(*b.get());
    }

    static bool equiv(const SafePtr<T>& b, const T* a) {
      return a->operator==(*b.get());
    }

    static bool equiv(const T* a, const T& b) {
      return a->operator==(b);
    }

#if !PTREQUIV_USE_TYPEID
    
    static bool equiv(const SafePtr<T>& a, const SafePtr<P>& b) {
      SafePtr<T> b_cast = dynamic_pointer_cast<T,P>(b);
      if (b_cast == 0)
        return false;
      else
        return a->operator==(*b_cast.get());
    }

    static bool equiv(const T* a, const SafePtr<P>& b) {
      SafePtr<T> b_cast = dynamic_pointer_cast<T,P>(b);
      if (b_cast == 0)
        return false;
      else
        return a->operator==(*b_cast.get());
    }

    static bool equiv(const T* a, const SafePtr<DGVertex>& b) {
      SafePtr<T> b_cast = dynamic_pointer_cast<T,DGVertex>(b);
      if (b_cast == 0)
        return false;
      else
        return a->operator==(*b_cast.get());
    }
    
#else

    static bool equiv(const SafePtr<T>& a, const SafePtr<P>& b) {
      if (a->typeid_ != b->typeid_)
        return false;
      else {
        SafePtr<T> b_cast = static_pointer_cast<T,P>(b);
        return a->operator==(*b_cast.get());
      }
    }

    static bool equiv(const T* a, const SafePtr<DGVertex>& b) {
      if (a->typeid_ != b->typeid_)
        return false;
      else {
#if PTREQUIV_USE_KEY_TO_COMPARE
  #if PTREQUIV_USE_INSTID
        return a->instid_ == b->instid_;
  #else
        return a->label() == b->label();
  #endif
#else
        SafePtr<T> b_cast = static_pointer_cast<T,DGVertex>(b);
        return a->operator==(*b_cast.get());
#endif
      }
    }

#endif

  };
    
    /*
     static bool equiv(const SafePtr<parent_type>& b, const SafePtr<T>& a) const {
       SafePtr<T> b_cast = dynamic_pointer_cast<T,parent_type>(b);
       if (b_cast == 0)
         return false;
       else
         return a->operator==(*b_cast.get())
     }
     */
        
};

#endif
