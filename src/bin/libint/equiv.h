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

#include <global_macros.h>
#include <smart_ptr.h>

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

  static bool equiv(const T& a, const T& b) { return a == b; }

  static bool equiv(const std::shared_ptr<T>& a, const std::shared_ptr<T>& b) {
    return a->operator==(*b.get());
  }

  static bool equiv(const T* a, const std::shared_ptr<T>& b) {
    return a->operator==(*b.get());
  }

  static bool equiv(const std::shared_ptr<T>& b, const T* a) {
    return a->operator==(*b.get());
  }

  static bool equiv(const T* a, const T& b) { return a->operator==(b); }

#if !PTREQUIV_USE_TYPEID

  static bool equiv(const std::shared_ptr<T>& a, const std::shared_ptr<P>& b) {
    std::shared_ptr<T> b_cast = std::dynamic_pointer_cast<T, P>(b);
    if (b_cast == 0)
      return false;
    else
      return a->operator==(*b_cast.get());
  }

  static bool equiv(const T* a, const std::shared_ptr<P>& b) {
    std::shared_ptr<T> b_cast = std::dynamic_pointer_cast<T, P>(b);
    if (b_cast == 0)
      return false;
    else
      return a->operator==(*b_cast.get());
  }

  static bool equiv(const T* a, const std::shared_ptr<DGVertex>& b) {
    std::shared_ptr<T> b_cast = std::dynamic_pointer_cast<T, DGVertex>(b);
    if (b_cast == 0)
      return false;
    else
      return a->operator==(*b_cast.get());
  }

#else

  static bool equiv(const std::shared_ptr<T>& a, const std::shared_ptr<P>& b) {
    if (a->typeid_ != b->typeid_)
      return false;
    else {
      std::shared_ptr<T> b_cast = std::static_pointer_cast<T, P>(b);
      return a->operator==(*b_cast.get());
    }
  }

  static bool equiv(const T* a, const std::shared_ptr<DGVertex>& b) {
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
      std::shared_ptr<T> b_cast = std::static_pointer_cast<T, DGVertex>(b);
      return a->operator==(*b_cast.get());
#endif
    }
  }

#endif
};

/*
 static bool equiv(const std::shared_ptr<parent_type>& b, const
 std::shared_ptr<T>& a) const { std::shared_ptr<T> b_cast =
 std::dynamic_pointer_cast<T,parent_type>(b); if (b_cast == 0) return false;
   else
     return a->operator==(*b_cast.get())
 }
 */

};  // namespace libint2

#endif
