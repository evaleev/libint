
#ifndef _libint2_src_bin_libint_singlstacktimpl_h_
#define _libint2_src_bin_libint_singlstacktimpl_h_

#define LOCAL_DEBUG 0

#include <singl_stack.h>
#include <purgeable.timpl.h>

namespace libint2 {

  template <class T, class KeyType>
    SingletonStack<T,KeyType>::SingletonStack(HashingFunction callback) :
    map_(), callback_(callback), next_instance_(0)
    {
      if (PurgingPolicy::purgeable()) { // if this stack contains objects that can be purged, add to the registry
        PurgeableStacks::Instance()->register_stack(this);
      }
    }

  template <class T, class KeyType>
    const typename SingletonStack<T,KeyType>::value_type&
    SingletonStack<T,KeyType>::find(const SafePtr<T>& obj)
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

  template <class T, class KeyType>
    const typename SingletonStack<T,KeyType>::value_type&
    SingletonStack<T,KeyType>::find(const key_type& key)
    {
      static value_type null_value(make_pair(InstanceID(0),SafePtr<T>()));
      typedef typename map_type::iterator miter;
      miter pos = map_.find(key);
      if (pos != map_.end()) {
        return (*pos).second;
      }
      else {
	return null_value;
      }
    }

  template <class T, class KeyType>
    void
    SingletonStack<T,KeyType>::remove(const SafePtr<T>& obj)
    {
      key_type key = ((obj.get())->*callback_)();

      typedef typename map_type::iterator miter;
      miter pos = map_.find(key);
      if (pos != map_.end()) {
        map_.erase(pos);
#if DEBUG || LOCAL_DEBUG
        std::cout << "Removed from stack " << obj->label() << std::endl;
#endif
      }
    }

  template <class T, class KeyType>
    void
    SingletonStack<T,KeyType>::purge()
    {
      for(iter_type i = map_.begin(); i!=map_.end();) {
        const T* v = i->second.second.get();
        if (PurgingPolicy::purge(v))
          // map::erase invalidates the iterator, increment it beforehand
          map_.erase(i++);
        else
          ++i;
      }
    }

};

#endif
