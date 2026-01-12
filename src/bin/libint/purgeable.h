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

#ifndef _libint2_src_bin_libint_purgeable_h_
#define _libint2_src_bin_libint_purgeable_h_

#include <dgvertex.h>

#include <boost/type_traits.hpp>
#include <vector>

namespace libint2 {

/** Determines whether an object should be purged from a stack. The default
   policy is to purge DGVector if it doesn't belong to any DirectedGraph.
  */
template <typename T>
struct DefaultPurgingPolicy {
  /// returns true if objects of this type can be purged
  static bool purgeable() {
    bool result = false;

    if (boost::is_base_of<DGVertex,
                          T>::value) {  // can only purge DGVertex objects
      result = true;
    }

    return result;
  }

  /// returns true if obj should be purged
  static bool purge(const T* ref) {
    bool result = false;

    try {
      const DGVertex* dgv_ptr = dynamic_cast<const DGVertex*>(ref);
      if (dgv_ptr->dg() == 0) result = true;
    } catch (...) {
    }

    return result;
  }
};

/**
 * PurgeableStack is a container that can be purged by calling purge() method.
 */
class AbstractPurgeableStack {
 public:
  virtual ~AbstractPurgeableStack() {}
  virtual void purge() = 0;
};

/**
 * PurgeableStack is an AbstractPurgeableStack that contains objects of type T.
 * Whether to purge an object is determined by calling Policy::purge(const T*)
 */
template <typename T, typename Policy = DefaultPurgingPolicy<T> >
class PurgeableStack : public AbstractPurgeableStack {
 protected:
  typedef Policy PurgingPolicy;

  virtual ~PurgeableStack() {}
};

/// Collection of AbstractPurgeableStack objects
class PurgeableStacks {
 public:
  typedef PurgeableStacks this_type;
  typedef AbstractPurgeableStack stack_type;

  static this_type* Instance();

  void purge();
  void register_stack(stack_type* stack);

 private:
  PurgeableStacks() {}
  static this_type* instance_;
  std::vector<stack_type*> stacks_;
};

};  // namespace libint2

#endif
