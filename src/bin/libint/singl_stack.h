
#ifndef _libint2_src_bin_libint_singlstack_h_
#define _libint2_src_bin_libint_singlstack_h_

#include <map>
#include <iostream>
#include <smart_ptr.h>
#include <key.h>
#include <hashable.h>
#include <purgeable.h>

namespace libint2 {

  class RecurrenceRelation;

  /**
     SingletonStack<T,KeyType> helps to implement Singleton-like objects of type T.
     SingletonStack maintains a map of keys of type KeyType to
     smart pointers to objects of type T. Keys are computed from T by
     calling a callback of type HashingFunction passed to the constructor to SingletonStack -- hence
     HashingFunction is a member function T which takes no arguments and returns a const KeyType&.
  */
  template <class T, class KeyType>
    class SingletonStack : public PurgeableStack<T>
    {
    public:
      typedef KeyType key_type;
      typedef SafePtr<T> data_type;
      /// Specifies type for the instance index variables
      typedef KeyTypes::InstanceID InstanceID;
      typedef std::pair<InstanceID,SafePtr<T> > value_type;
      typedef std::map<key_type,value_type> map_type;
      /// use iter_type objects to iterate over the stack
      typedef typename map_type::iterator iter_type;
      /// const version of iter_type
      typedef typename map_type::const_iterator citer_type;
      /// hashing function returns keys as key_return_type
      typedef typename KeyTraits<key_type>::ReturnType key_return_type;
      /// Specifies the type of callback which computes hashes
      typedef key_return_type (T::* HashingFunction)() const;
      /// PurgingPolicy determines whether and which objects on this stack are obsolete and can be removed
      typedef typename PurgeableStack<T>::PurgingPolicy PurgingPolicy;


      /// callback to compute hash values is the only parameter
      SingletonStack(HashingFunction callback);
      virtual ~SingletonStack() {}

      /** Returns the pointer to the unique instance of object obj.
          find() computes obj->*callback_(), searches it in hstack_,
          if found -- returns the pointer to the corresponding object on ostack_,
          otherwise pushes obj to the end of ostack_ and returns obj.
      */
      const value_type& find(const SafePtr<T>& obj);
      /** Returns the pointer to the unique instance of object corresponding to key.
          if found returns the pointer to the corresponding object on ostack_,
          else returns a null value_type.
      */
      const value_type& find(const key_type& key);

      /** Searches for obj on the stack and, if found, removes the unique instance
      */
      void remove(const SafePtr<T>& obj);

      /** Returns iterator to the beginning of the stack */
      citer_type begin() const { return map_.begin(); }
      /** Returns iterator to the end of the stack */
      citer_type end() const { return map_.end(); }

      // Implementation of PurgeableStack::purge()
      void purge();

    private:
      map_type map_;
      HashingFunction callback_;
      InstanceID next_instance_;
    };

};

#endif
