/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint compiler.
 *
 *  Libint compiler is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint compiler is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint compiler.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_bin_libint_singlstack_h_
#define _libint2_src_bin_libint_singlstack_h_

#include <hashable.h>
#include <key.h>
#include <purgeable.h>
#include <smart_ptr.h>

#include <iostream>
#include <map>

namespace libint2 {

class RecurrenceRelation;

/**
   SingletonStack<T,KeyType> helps to implement Singleton-like objects of type
   T. SingletonStack maintains a map of keys of type KeyType to smart pointers
   to objects of type T. Keys are computed from T by calling a callback of type
   HashingFunction passed to the constructor to SingletonStack -- hence
   HashingFunction is a member function T which takes no arguments and returns a
   const KeyType&.
 */
template <class T, class KeyType>
class SingletonStack : public PurgeableStack<T> {
 public:
  typedef KeyType key_type;
  typedef std::shared_ptr<T> data_type;
  /// Specifies type for the instance index variables
  typedef KeyTypes::InstanceID InstanceID;
  typedef std::pair<InstanceID, std::shared_ptr<T> > value_type;
  typedef std::map<key_type, value_type> map_type;
  /// use iter_type objects to iterate over the stack
  typedef typename map_type::iterator iter_type;
  /// const version of iter_type
  typedef typename map_type::const_iterator citer_type;
  /// hashing function returns keys as key_return_type
  typedef typename KeyTraits<key_type>::ReturnType key_return_type;
  /// Specifies the type of callback which computes hashes
  typedef key_return_type (T::*HashingFunction)() const;
  /// PurgingPolicy determines whether and which objects on this stack are
  /// obsolete and can be removed
  typedef typename PurgeableStack<T>::PurgingPolicy PurgingPolicy;

  /// callback to compute hash values is the only parameter
  SingletonStack(HashingFunction callback)
      : map_(), callback_(callback), next_instance_(0) {
    if (PurgingPolicy::purgeable()) {  // if this stack contains objects that
                                       // can be purged, add to the registry
      PurgeableStacks::Instance()->register_stack(this);
    }
  }

  virtual ~SingletonStack() {}

  /** Returns the pointer to the unique instance of object obj.
        find() computes obj->*callback_(), searches it in hstack_,
        if found -- returns the pointer to the corresponding object on ostack_,
        otherwise pushes obj to the end of ostack_ and returns obj.
   */
  const value_type& find(const std::shared_ptr<T>& obj) {
    key_type key = ((obj.get())->*callback_)();

    typedef typename map_type::iterator miter;
    miter pos = map_.find(key);
    if (pos != map_.end()) {
#if DEBUG || LOCAL_DEBUG
      std::cout << "SingletonStack::find -- " << obj->label()
                << " already found" << std::endl;
#endif
      return (*pos).second;
    } else {
      value_type result(next_instance_++, obj);
      map_[key] = result;
#if DEBUG || LOCAL_DEBUG
      std::cout << "SingletonStack::find -- " << obj->label()
                << " is new (instid_ = " << next_instance_ - 1 << ")"
                << std::endl;
#endif
      return map_[key];
    }
  }

  /** Returns the pointer to the unique instance of object corresponding to key.
        if found returns the pointer to the corresponding object on ostack_,
        else returns a null value_type.
   */
  const value_type& find(const key_type& key) {
    static value_type null_value(
        make_pair(InstanceID(0), std::shared_ptr<T>()));
    typedef typename map_type::iterator miter;
    miter pos = map_.find(key);
    if (pos != map_.end()) {
      return (*pos).second;
    } else {
      return null_value;
    }
  }

  /** Returns the pointer to the unique instance of object corresponding to the
     hashed_key. if found returns the pointer to the corresponding object on
     ostack_, else returns a null value_type.
   */
  const value_type& find_hashed(const InstanceID& hashed_key) const {
    static value_type null_value(
        make_pair(InstanceID(0), std::shared_ptr<T>()));
    for (auto& i : map_) {
      if (i.second.first == hashed_key) return i.second;
    }
    return null_value;
  }

  /** Searches for obj on the stack and, if found, removes the unique instance
   */
  void remove(const std::shared_ptr<T>& obj) {
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

  /** Returns iterator to the beginning of the stack */
  citer_type begin() const { return map_.begin(); }
  /** Returns iterator to the end of the stack */
  citer_type end() const { return map_.end(); }

  // Implementation of PurgeableStack::purge()
  void purge() override {
    for (iter_type i = map_.begin(); i != map_.end();) {
      const T* v = i->second.second.get();
      if (PurgingPolicy::purge(v))
        // map::erase invalidates the iterator, increment it beforehand
        map_.erase(i++);
      else
        ++i;
    }
  }

 private:
  map_type map_;
  HashingFunction callback_;
  InstanceID next_instance_;
};

};  // namespace libint2

#endif
