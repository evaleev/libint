/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint library.
 *
 *  Libint library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_lib_libint_util_compressedpair_h_
#define _libint2_src_lib_libint_util_compressedpair_h_

namespace libint2 {
namespace detail {

// a bare-bones implementation that assumes T1 and T2 are classes, and T1 != T2
template <class T1, class T2>
class compressed_pair : private T1, private T2 {
 public:
  typedef T1 first_type;
  typedef T2 second_type;
  typedef typename std::add_const<first_type>::type first_const_type;
  typedef typename std::add_const<second_type>::type second_const_type;
  typedef typename std::add_lvalue_reference<first_type>::type first_reference;
  typedef
      typename std::add_lvalue_reference<second_type>::type second_reference;
  typedef typename std::add_lvalue_reference<first_const_type>::type
      first_const_reference;
  typedef typename std::add_lvalue_reference<second_const_type>::type
      second_const_reference;
  typedef typename std::add_rvalue_reference<first_type>::type
      first_rvalue_reference;
  typedef typename std::add_rvalue_reference<second_type>::type
      second_rvalue_reference;

  compressed_pair() = default;
  compressed_pair(const first_type& x, const second_type& y)
      : first_type(x), second_type(y) {}
  explicit compressed_pair(const first_type& x) : first_type(x) {}
  explicit compressed_pair(const second_type& y) : second_type(y) {}

  compressed_pair(const compressed_pair& other) = default;
  compressed_pair(compressed_pair&& other)
      : first_type(std::move(other.first_rvalref())),
        second_type(std::move(other.second_rvalref())) {}
  compressed_pair& operator=(const compressed_pair&) = default;
  compressed_pair& operator=(compressed_pair&& other) {
    this->first() = std::move(other.first_rvalref());
    this->second() = std::move(other.second_rvalref());
    return *this;
  }

  first_reference first() { return static_cast<first_reference>(*this); }
  first_const_reference first() const {
    return static_cast<first_const_reference>(*this);
  }
  first_rvalue_reference first_rvalref() {
    return static_cast<first_rvalue_reference>(*this);
  }

  second_reference second() { return static_cast<second_reference>(*this); }
  second_const_reference second() const {
    return static_cast<second_const_reference>(*this);
  }
  second_rvalue_reference second_rvalref() {
    return static_cast<second_rvalue_reference>(*this);
  }

  void swap(compressed_pair& other) {
    swap(this->first(), other.first());
    swap(this->second(), other.second());
  }
};

template <class T1, class T2>
compressed_pair<T1, T2> make_compressed_pair(const T1& x, const T2& y) {
  return compressed_pair<T1, T2>(x, y);
}

}  // namespace detail
}  // namespace libint2
#endif  // header guard
