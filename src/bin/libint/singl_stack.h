
#include <map>
#include <iostream>
#include <smart_ptr.h>

#define LOCAL_DEBUG 0

#ifndef _libint2_src_bin_libint_singlstack_h_
#define _libint2_src_bin_libint_singlstack_h_

namespace libint2 {

  class RecurrenceRelation;
  
  /// This is used to maintain some information about Generalized Singletons
  struct GSingletonTrait {
    typedef unsigned long int InstanceID;
  };
  
  /**
     SingletonStack<T,HashType,P> helps to implement Singleton-like objects of type T.
     SingletonStack maintains a map of keys of type HashType to
     smart pointers to objects of type T. Keys are computed from T by
     calling a callback of type HashCallback passed to the constructor to SingletonStack -- hence
     HashCallback is a member function T which takes no arguments and returns a const HashType&.
  */
  template <class T, class HashType>
    class SingletonStack
    {
    public:
      typedef HashType key_type;
      typedef SafePtr<T> data_type;
      /// Specifies type for the instance index variables
      typedef GSingletonTrait::InstanceID InstanceID;
      typedef std::pair<InstanceID,SafePtr<T> > value_type;
      typedef std::map<key_type,value_type> map_type;
      /// use iter_type objects to iterate over the stack
      typedef typename map_type::iterator iter_type;
      /// const version of iter_type
      typedef typename map_type::const_iterator citer_type;
      /// Specifies the type of callback which computes hashes
      typedef const key_type& (T::* HashCallback)() const;

      /// callback to compute hash values is the only parameter
      SingletonStack(HashCallback callback);
      ~SingletonStack() {}

      /** Returns the pointer to the unique instance of object obj.
          find() computes obj->*callback_(), searches it in hstack_,
          if found -- returns the pointer to the corresponding object on ostack_,
          otherwise pushes obj to the end of ostack_ and returns obj.
      */
      const value_type& find(const SafePtr<T>& obj);
      
      /** Returns iterator to the beginning of the stack */
      citer_type begin() const { return map_.begin(); }
      /** Returns iterator to the end of the stack */
      citer_type end() const { return map_.end(); }

    private:
      map_type map_;
      HashCallback callback_;
      InstanceID next_instance_;
    };

  template <class T, class HashType>
    SingletonStack<T,HashType>::SingletonStack(HashCallback callback) :
    map_(), callback_(callback), next_instance_(0)
    {
    }

  template <class T, class HashType>
    const typename SingletonStack<T,HashType>::value_type&
    SingletonStack<T,HashType>::find(const SafePtr<T>& obj)
    {
      key_type key = ((obj.get())->*callback_)();

      typedef typename map_type::iterator miter;
      miter pos = map_.find(key);
      if (pos != map_.end()) {
#if DEBUG || LOCAL_DEBUG
        std::cout << "SingletonStack::find -- " << obj->label() << " already found" << std::endl;
#endif
        return (*pos).second;
      }
      else {
        value_type result(next_instance_++,obj);
        map_[key] = result;
#if DEBUG || LOCAL_DEBUG
        std::cout << "SingletonStack::find -- " << obj->label() << " is new (instid_ = " << next_instance_-1 << ")" << std::endl;
#endif
        return map_[key];
      }
    }

};

#endif
