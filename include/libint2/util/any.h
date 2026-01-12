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

#ifndef _libint2_include_libint2_util_any_h_
#define _libint2_include_libint2_util_any_h_

#include <cassert>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <utility>

// Include C++17 any header, if available AND functional
#if __cplusplus >= 201703L
// macos < 10.14 do not have any_cast in their libc++
#include <ciso646>  // see https://stackoverflow.com/a/31658120
#if defined(_LIBCPP_VERSION) && defined(__APPLE__)
#include <Availability.h>
#ifdef __MAC_OS_X_VERSION_MIN_ALLOWED
#if __MAC_OS_X_VERSION_MIN_ALLOWED >= 10140
#define LIBINT_HAS_CXX17_ANY
#endif  //  10.14 or later
#endif  // have macos version
#else   // libc++ on macos
#define LIBINT_HAS_CXX17_ANY
#endif  // libc++ on macos
#endif  // c++17

#ifdef LIBINT_HAS_CXX17_ANY
#include <any>
#endif

namespace libint2 {

// prefer std::any, if available
#ifdef LIBINT_HAS_CXX17_ANY
using std::any;
using std::any_cast;
using std::bad_any_cast;
#else

namespace detail {
// true if decayed T is Base, or is derived from it
template <typename Base, typename T>
using disable_if_same_or_derived = typename std::enable_if<
    !std::is_base_of<Base, typename std::decay<T>::type>::value>::type;
};  // namespace detail

/// a partial C++17 std::any implementation (and less efficient than can be)
class any {
 public:
  // this is constexpr in the standard
  any() : impl_(nullptr) {}
  any(const any& other) : impl_(other.impl_->clone()) {}
  any(any&& other) = default;
  template <typename ValueType,
            typename = detail::disable_if_same_or_derived<any, ValueType> >
  any(ValueType&& value)
      : impl_(new impl<typename std::decay<ValueType>::type>(
            std::forward<ValueType>(value))) {}
  ~any() = default;

  any& operator=(const any& rhs) {
    impl_ = decltype(impl_)(rhs.impl_->clone());
    return *this;
  }
  any& operator=(any&& rhs) {
    impl_ = std::move(rhs.impl_);
    return *this;
  }
  template <typename ValueType,
            typename = detail::disable_if_same_or_derived<any, ValueType> >
  any& operator=(ValueType&& rhs) {
    impl_ = decltype(impl_)(new impl<typename std::decay<ValueType>::type>(
        std::forward<ValueType>(rhs)));
    return *this;
  }

  template <class ValueType, class... Args>
  typename std::decay<ValueType>::type& emplace(Args&&... args) {
    reset();
    impl_ = new impl<typename std::decay<ValueType>::type>(
        std::forward<Args>(args)...);
    return (impl_->cast_static<typename std::decay<ValueType>::type>()->value);
  }
  template <class ValueType, class U, class... Args>
  typename std::decay<ValueType>::type& emplace(std::initializer_list<U> il,
                                                Args&&... args) {
    reset();
    impl_ = new impl<typename std::decay<ValueType>::type>(
        il, std::forward<Args>(args)...);
    return (impl_->cast_static<typename std::decay<ValueType>::type>()->value);
  }

  void reset() { impl_.reset(); }

  void swap(any& other) { std::swap(impl_, other.impl_); }

  bool has_value() const { return static_cast<bool>(impl_); }

  const std::type_info& type() const {
    if (has_value())
      return impl_->type();
    else
      return typeid(void);
  }

 private:
  template <typename T>
  struct impl;

  struct impl_base {
    virtual ~impl_base() {}
    virtual impl_base* clone() const = 0;

    virtual const std::type_info& type() const = 0;

    // static if NDEBUG is defined, dynamic otherwise
    template <typename T>
    impl<T>* cast() {
#ifndef NDEBUG
      return this->cast_static<T>();
#else
      return dynamic_cast<impl<T>*>(this);
#endif
    }
    // static always
    template <typename T>
    impl<T>* cast_static() {
      return static_cast<impl<T>*>(this);
    }
  };
  template <typename T>
  struct impl : public impl_base {
    template <typename U>
    explicit impl(U&& v) : value(std::forward<U>(v)) {}
    impl_base* clone() const override { return new impl{value}; }

    const std::type_info& type() const override { return typeid(T); }

    T value;
  };

  template <typename ValueType>
  friend typename std::decay<ValueType>::type* any_cast(any* operand);
  template <typename ValueType>
  friend const typename std::decay<ValueType>::type* any_cast(
      const any* operand);

  template <typename ValueType>
  typename std::decay<ValueType>::type* value_ptr() {
    return &(impl_->cast_static<typename std::decay<ValueType>::type>()->value);
  }

  template <typename ValueType>
  const typename std::decay<ValueType>::type* value_ptr() const {
    return &(impl_->cast_static<typename std::decay<ValueType>::type>()->value);
  }

  std::unique_ptr<impl_base> impl_;
};

class bad_any_cast : public std::bad_cast {
 public:
  bad_any_cast() = default;
  virtual ~bad_any_cast() {}
  virtual const char* what() const noexcept { return "Bad any_cast"; }
};

template <typename ValueType>
typename std::decay<ValueType>::type* any_cast(any* operand) {
  if (operand->type() == typeid(typename std::decay<ValueType>::type))
    return operand->value_ptr<typename std::decay<ValueType>::type>();
  else
    return nullptr;
}

template <typename ValueType>
const typename std::decay<ValueType>::type* any_cast(const any* operand) {
  if (operand->type() == typeid(typename std::decay<ValueType>::type))
    return operand->value_ptr<typename std::decay<ValueType>::type>();
  else
    return nullptr;
}

template <typename ValueType>
ValueType any_cast(const any& operand) {
  const auto* cast_ptr =
      any_cast<typename std::decay<ValueType>::type>(&operand);
  if (cast_ptr != nullptr) return *cast_ptr;
  throw bad_any_cast();
}

template <typename ValueType>
ValueType any_cast(any& operand) {
  auto* cast_ptr = any_cast<typename std::decay<ValueType>::type>(&operand);
  if (cast_ptr != nullptr) return *cast_ptr;
  throw bad_any_cast();
}
#endif  // C++17

}  // namespace libint2

#endif  // header guard
