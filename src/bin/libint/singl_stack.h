
#include <vector>
#include <algorithm>
#include <smart_ptr.h>

#ifndef _libint2_src_bin_libint_singlstack_h_
#define _libint2_src_bin_libint_singlstack_h_

namespace libint2 {

  /**
     SingletonStack<T,HashType> helps to implement Singleton-like objects of type T.
     SingletonStack maintains a container of smart pointers to objects of type T and a
     corresponding container of hashes of type HashType. Hashes are computed from T by
     calling a callback of type HashCallback passed to the constructor to SingletonStack -- hence
     HashCallback is a member function of T which takes on arguments and returns HashType.
  */
  template <class T, class HashType>
    class SingletonStack
    {
    public:
      typedef std::vector< SafePtr<T> > ostack_type;
      typedef std::vector< HashType > hstack_type;
      /// Specifies the type of callback which computes hashes
      typedef HashType (T::* HashCallback)() const;

      /// callback to compute hash values is the only parameter
      SingletonStack(HashCallback callback);
      ~SingletonStack() {}

      /** Returns the pointer to the unique instance of object obj.
          find() computes obj->*callback_(), searches it in hstack_,
          if found -- returns the pointer to the corresponding object on ostack_,
          otherwise pushes obj to the end of ostack_ and returns obj.
      */
      SafePtr<T> find(const SafePtr<T>& obj);

    private:
      ostack_type ostack_;
      hstack_type hstack_;
      HashCallback callback_;
    };

  template <class T, class HashType>
    SingletonStack<T,HashType>::SingletonStack(HashCallback callback) :
    ostack_(), hstack_(), callback_(callback)
    {
    }

  template <class T, class HashType>
    SafePtr<T>
    SingletonStack<T,HashType>::find(const SafePtr<T>& obj)
    {
      HashType hash = ((obj.get())->*callback_)();
      typedef typename hstack_type::iterator hiter;
      hiter begin = hstack_.begin();
      hiter end = hstack_.end();
      hiter pos = std::find(begin,end,hash);
      if (pos != end) {
        return ostack_[pos-begin];
      }

      ostack_.push_back(obj);
      hstack_.push_back(hash);
      return obj;
    }

};

#endif
