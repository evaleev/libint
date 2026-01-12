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

#ifndef _libint2_include_libint2_util_optional_h_
#define _libint2_include_libint2_util_optional_h_

#if __cplusplus < 201703L

#include <initializer_list>  // For initializer_list
#include <new>               // For placement new
#include <type_traits>       // For is_trivially_copyable, etc.
#include <utility>           // For std::move

#else

#include <optional>

#endif

namespace libint2 {

#if __cplusplus < 201703L

struct nullopt_t {
  constexpr explicit nullopt_t(int) {}
};
constexpr nullopt_t nullopt{0};

struct bad_optional_access : public std::exception {
  bad_optional_access() noexcept = default;
  bad_optional_access(const bad_optional_access &other) noexcept = default;
  bad_optional_access &operator=(const bad_optional_access &other) noexcept =
      default;
};

struct in_place_t {
  explicit in_place_t() = default;
};

template <typename T>
class optional {
 private:
  bool has_value_;
  typename std::aligned_storage<sizeof(T), alignof(T)>::type storage_;

 public:
  // Constructors
  constexpr optional() noexcept : has_value_(false) {}

  constexpr optional(nullopt_t) noexcept : has_value_(false) {}

  optional(const optional &other) : has_value_(other.has_value_) {
    if (has_value_) {
      new (&storage_) T(*other);
    }
  }

  optional(optional &&other) noexcept : has_value_(other.has_value_) {
    if (has_value_) {
      new (&storage_) T(std::move(*other));
    }
    other.has_value_ = false;
  }

  template <typename... Args>
  optional(in_place_t, Args &&...args) : has_value_(true) {
    new (&storage_) T(std::forward<Args>(args)...);
  }

  template <typename U, typename... Args>
  optional(in_place_t, std::initializer_list<U> il, Args &&...args)
      : has_value_(true) {
    new (&storage_) T(il, std::forward<Args>(args)...);
  }

  template <
      typename U = T,
      typename std::enable_if<
          std::is_convertible<U, T>::value &&
              !std::is_same<optional, typename std::decay<U>::type>::value,
          int>::type = 0>
  optional(U &&value) : has_value_(true) {
    new (&storage_) T(std::forward<U>(value));
  }

  // Destructor
  ~optional() {
    if (has_value_) {
      destroy();
    }
  }

  // Assignment operators
  optional &operator=(const optional &other) {
    if (this != &other) {
      if (has_value_ && !other.has_value_) {
        destroy();
        has_value_ = false;
      } else if (!has_value_ && other.has_value_) {
        new (&storage_) T(*other);
        has_value_ = true;
      } else if (has_value_ && other.has_value_) {
        *this = other;
      }
    }
    return *this;
  }

  optional &operator=(optional &&other) noexcept {
    if (this != &other) {
      if (has_value_) {
        destroy();
      }
      if (other.has_value_) {
        new (&storage_) T(std::move(*other));
        has_value_ = true;
      } else {
        has_value_ = false;
      }
      other.has_value_ = false;
    }
    return *this;
  }

  template <
      typename U = T,
      typename std::enable_if<std::is_convertible<U, T>::value, int>::type = 0>
  optional &operator=(U &&value) {
    if (has_value_) {
      destroy();
    }
    new (&storage_) T(std::forward<U>(value));
    has_value_ = true;
    return *this;
  }

  // Observers
  bool has_value() const noexcept { return has_value_; }

  T &operator*() {
    if (!has_value_) throw bad_optional_access();
    return *get_ptr();
  }

  const T &operator*() const {
    if (!has_value_) throw bad_optional_access();
    return *get_ptr();
  }

  T *operator->() {
    if (!has_value_) throw bad_optional_access();
    return get_ptr();
  }

  const T *operator->() const {
    if (!has_value_) throw bad_optional_access();
    return get_ptr();
  }

  T &value() {
    if (!has_value_) throw bad_optional_access();
    return *get_ptr();
  }

  const T &value() const {
    if (!has_value_) throw bad_optional_access();
    return *get_ptr();
  }

  template <typename U>
  T value_or(U &&default_value) const & {
    return has_value_ ? *get_ptr()
                      : static_cast<T>(std::forward<U>(default_value));
  }

  template <typename U>
  T value_or(U &&default_value) && {
    return has_value_ ? std::move(*get_ptr())
                      : static_cast<T>(std::forward<U>(default_value));
  }

  // Modifiers
  void reset() noexcept {
    if (has_value_) {
      destroy();
      has_value_ = false;
    }
  }

  template <typename... Args>
  T &emplace(Args &&...args) {
    reset();
    new (&storage_) T(std::forward<Args>(args)...);
    has_value_ = true;
    return *get_ptr();
  }

  template <typename U, typename... Args>
  T &emplace(std::initializer_list<U> il, Args &&...args) {
    reset();
    new (&storage_) T(il, std::forward<Args>(args)...);
    has_value_ = true;
    return *get_ptr();
  }

  void swap(optional &other) noexcept {
    if (this != &other) {
      std::swap(has_value_, other.has_value_);
      if (has_value_ && other.has_value_) {
        std::swap(*get_ptr(), *other.get_ptr());
      } else if (has_value_ && !other.has_value_) {
        new (&other.storage_) T(std::move(*this));
        destroy();
        std::swap(storage_, other.storage_);
      } else if (!has_value_ && other.has_value_) {
        new (&storage_) T(std::move(*other));
        destroy();
        std::swap(storage_, other.storage_);
      }
    }
  }

 private:
  T *get_ptr() { return reinterpret_cast<T *>(&storage_); }
  const T *get_ptr() const { return reinterpret_cast<const T *>(&storage_); }

  void destroy() { get_ptr()->~T(); }
};

template <typename T>
void swap(optional<T> &a, optional<T> &b) noexcept {
  a.swap(b);
}

template <typename T>
bool operator==(const optional<T> &x, const optional<T> &y) {
  return (!x && !y) || (x && y && *x == *y);
}

template <typename T>
bool operator!=(const optional<T> &x, const optional<T> &y) {
  return !(x == y);
}

template <typename T>
bool operator<(const optional<T> &x, const optional<T> &y) {
  return y && (!x || *x < *y);
}

template <typename T>
bool operator>(const optional<T> &x, const optional<T> &y) {
  return y < x;
}

template <typename T>
bool operator<=(const optional<T> &x, const optional<T> &y) {
  return !(x > y);
}

template <typename T>
bool operator>=(const optional<T> &x, const optional<T> &y) {
  return !(x < y);
}

#else

using std::nullopt;
using std::optional;

#endif

}  // namespace libint2

#endif  // _libint2_include_libint2_util_optional_h_
