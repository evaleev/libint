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

#include <bfset.h>
#include <global_macros.h>
#include <smart_ptr.h>

#ifndef _libint2_src_bin_libint_traits_h_
#define _libint2_src_bin_libint_traits_h_

namespace libint2 {

template <typename T>
struct StorageTraits {
  typedef std::shared_ptr<T> StorageType;
  enum { StoredAsPtr = true };
  static const T& const_ref(const StorageType& s) { return *s; };
};

#if USE_BRAKET_H
template <>
struct StorageTraits<CGShell> {
  typedef CGShell StorageType;
  enum { StoredAsPtr = false };
  static const CGShell& const_ref(const StorageType& s) { return s; };
};

template <>
struct StorageTraits<CGF> {
  typedef CGF StorageType;
  enum { StoredAsPtr = false };
  static const CGF& const_ref(const StorageType& s) { return s; };
};

template <CartesianAxis Axis>
struct StorageTraits<CGShell1d<Axis> > {
  typedef CGShell1d<Axis> StorageType;
  enum { StoredAsPtr = false };
  static const CGShell1d<Axis>& const_ref(const StorageType& s) { return s; };
};

template <CartesianAxis Axis>
struct StorageTraits<CGF1d<Axis> > {
  typedef CGF1d<Axis> StorageType;
  enum { StoredAsPtr = false };
  static const CGF1d<Axis>& const_ref(const StorageType& s) { return s; };
};
#endif

///////////

/// Converts Base to a type of the same signature as Ref. For example, if Ref is
/// std::shared_ptr<T> then Base is converted to std::shared_ptr<Base>
template <typename Ref, typename Base>
struct ReturnTypeAnalog {
  typedef const Base& result;
};
template <typename Ref, typename Base>
struct ReturnTypeAnalog<std::shared_ptr<Ref>, Base> {
  typedef std::shared_ptr<Base> result;
};

///////////

template <typename T>
struct TypeTraits {
  /// By default, use std::shared_ptr to manage these objects
  typedef typename StorageTraits<T>::StorageType StorageType;
  /// Whether stored as a pointer
  enum { StoredAsPtr = StorageTraits<T>::StoredAsPtr };
  /// Convert an object of StorageType to const T&
  static const T& const_ref(const StorageType& s) {
    return StorageTraits<T>::const_ref(s);
  }
};
};  // namespace libint2

#endif
