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

#ifndef _libint2_src_bin_libint_policy_h_
#define _libint2_src_bin_libint_policy_h_

#include <smart_ptr.h>
#include <traits.h>

#include <boost/type_traits/is_same.hpp>
#include <vector>

namespace libint2 {

/**
  ExistsDefaultSubobjAllocator<T,cond> is a helper class that defines the
  default subobj allocator for T. Such allocator is defined only if cond is
  true.
*/
template <class T, bool exists>
struct ExistsDefaultSubobjAllocator;

template <class T>
struct ExistsDefaultSubobjAllocator<T, true> {
  typedef T obj_type;
  typedef typename TypeTraits<obj_type>::StorageType obj_stype;
  typedef typename TypeTraits<obj_type>::StorageType subobj_stype;
  typedef typename obj_type::iter_type subobj_type;

  static void default_init_subobj(const obj_stype& obj,
                                  std::vector<subobj_stype>& subobj) {
    subobj.push_back(obj);
  }
  static void default_dealloc_subobj(std::vector<subobj_stype>& subobj) {}
};

/**
    StdLibintTDPolicy<T> is the default type-specific policy.
*/
template <class T>
struct StdLibintTDPolicy {
  typedef T obj_type;
  typedef typename obj_type::iter_type subobj_type;
  /// how these objects are stored
  typedef typename TypeTraits<obj_type>::StorageType obj_stype;
  /// how these subobjects are stored
  typedef typename TypeTraits<subobj_type>::StorageType subobj_stype;

  /// This function allocates subobj of obj (e.g. basis functions contained in a
  /// shell)
  static void init_subobj(const obj_stype& obj,
                          std::vector<subobj_stype>& subobj) {
    // If types are not the same then this function should not be used -- user
    // must provide a specialization
    ExistsDefaultSubobjAllocator<T,
                                 boost::is_same<obj_type, subobj_type>::value>::
        default_init_subobj(obj, subobj);
  }
  static void dealloc_subobj(std::vector<subobj_stype>& subobj) {
    // If types are not the same then this function should not be used -- user
    // must provide a specialization
    ExistsDefaultSubobjAllocator<T, boost::is_same<obj_type, subobj_type>::
                                        value>::default_dealloc_subobj(subobj);
  }
};

/**
    StdLibintTIPolicy is the default type-independent policy.
 */
class StdLibintTIPolicy {
  /// This controls whether can unroll an integral set
  static unsigned int max_set_size_to_unroll_;

 public:
  StdLibintTIPolicy() {}

  virtual void set_max_set_size_to_unroll(unsigned int i) {
    max_set_size_to_unroll_ = i;
  }

  virtual unsigned int max_set_size_to_unroll() const {
    return max_set_size_to_unroll_;
  }
};

/**
  Policy<T, TDPol, TIPol> defines a policy for type T as a combination of
 type-independent (TIPol) policies and type-specific (TDPol) policies.
 */
#if CXX_ALLOWS_DEFPARAMTEMPLATE_AS_TEMPTEMPPARAM
template <class T, class TIPol = StdLibintTIPolicy,
          template <class> class TDPol = StdLibintTDPolicy>
class Policy : public TDPol<T>, public TIPol {
#else
#define TDPol StdLibintTDPolicy
#define TIPol StdLibintTIPolicy
template <class T>
class Policy : public TDPol<T>, public TIPol {
#endif
 public:
  /// how these objects are stored
  typedef typename TDPol<T>::obj_stype obj_stype;
  /// how these subobjects are stored
  typedef typename TDPol<T>::subobj_stype subobj_stype;

  /*
    typedef typename TDPol<T>::obj_type obj_type;
    typedef typename obj_type::iter_type subobj_type;

    static void init_subobj(const std::shared_ptr<obj_type>& obj, const vector<
    std::shared_ptr<subobj_type> >& subobj)
    {
      TDPol<T>::init_subobj(obj,subobj);
    }

    static void dealloc_subobj(const vector< std::shared_ptr<subobj_type> >&
    subobj)
    {
      TDPol<T>::dealloc_subobj(subobj);
    }
   */

 private:
  /// Returns whether to unroll an IntegralSet
  bool can_unroll_intset(const std::shared_ptr<T>& iset) {
    return iset->set_size() <= TIPol::max_set_size_to_unroll();
  }
};

};  // namespace libint2

#endif
