
#include <map>
#include <smart_ptr.h>

#ifndef _libint2_src_bin_libint_singlstack_h_
#define _libint2_src_bin_libint_singlstack_h_

namespace libint2 {

  /**
     SingletonStack<T,HashType> helps to implement Singleton-like objects of type T.
     SingletonStack maintains a map of keys of type HashType to
     smart pointers to objects of type T. Keys are computed from T by
     calling a callback of type HashCallback passed to the constructor to SingletonStack -- hence
     HashCallback is a member function of T which takes on arguments and returns HashType.
  */
  template <class T, class HashType>
    class SingletonStack
    {
    public:
      typedef HashType key_type;
      typedef SafePtr<T> data_type;
      typedef std::map<key_type,data_type> map_type;
      /// Specifies the type of callback which computes hashes
      typedef key_type (T::* HashCallback)() const;

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
      map_type map_;
      HashCallback callback_;
    };

  template <class T, class HashType>
    SingletonStack<T,HashType>::SingletonStack(HashCallback callback) :
    map_(), callback_(callback)
    {
    }

  template <class T, class HashType>
    SafePtr<T>
    SingletonStack<T,HashType>::find(const SafePtr<T>& obj)
    {
      key_type key = ((obj.get())->*callback_)();

      typedef typename map_type::iterator miter;
      miter pos = map_.find(key);
      if (pos != map_.end()) {
        return (*pos).second;
      }
      else {
        map_[key] = obj;
        return obj;
      }
    }

};

#endif
