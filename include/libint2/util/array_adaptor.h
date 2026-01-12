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

#ifndef INCLUDE_LIBINT2_UTIL_ARRAY_ADAPTOR_H_
#define INCLUDE_LIBINT2_UTIL_ARRAY_ADAPTOR_H_

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <type_traits>
#include <vector>

namespace libint2 {
namespace detail {

/// allocator that uses an externally-managed stack-allocated array for
/// allocations up to max_size, for larger allocations uses heap.
template <class T, std::size_t N>
class ext_stack_allocator {
 public:
  using value_type = T;
  using pointer = T*;
  using difference_type =
      typename std::pointer_traits<pointer>::difference_type;
  using size_type = typename std::make_unsigned<difference_type>::type;
  using propagate_on_container_copy_assignment = std::true_type;
  using propagate_on_container_move_assignment = std::true_type;

  static auto constexpr size = N;
  typedef T array_type[N];

 private:
  T* stack_;  // stack-allocated storage
  T* free_;   // ptr to the first free element on stack

 public:
  ext_stack_allocator() noexcept : stack_(nullptr), free_(stack_) {}
  ext_stack_allocator(const ext_stack_allocator& other) = default;
  ext_stack_allocator(ext_stack_allocator&& other) noexcept
      : stack_(other.stack_), free_(other.free_) {}
  ext_stack_allocator& operator=(const ext_stack_allocator& other) = default;
  ext_stack_allocator& operator=(ext_stack_allocator&& other) noexcept {
    stack_ = other.stack_;
    free_ = other.free_;
    return *this;
  }

  explicit ext_stack_allocator(array_type& array) noexcept
      : stack_(&array[0]), free_(stack_) {}
  template <typename U,
            typename = typename std::enable_if<std::is_same<const U, T>::value>>
  explicit ext_stack_allocator(U (&array)[N]) noexcept
      : stack_(const_cast<T*>(&array[0])), free_(stack_) {}

  template <class _Up>
  struct rebind {
    using other = ext_stack_allocator<_Up, N>;
  };

  T* allocate(std::size_t n) {
    assert(stack_ != nullptr && "array_view_allocator not initialized");
    if (stack_ + N - free_ >=
        static_cast<std::ptrdiff_t>(n)) {  // have free space on stack
      const auto result = free_;
      free_ += n;
      return result;
    } else {
      return new T[n];
    }
  }
  void deallocate(T* p, std::size_t n) noexcept {
    if (pointer_on_stack(p)) {
      assert(p + n == free_ && "stack deallocation out of order");
      free_ -= n;
    } else {
      delete[] p;
    }
  }

  template <class T1, std::size_t N1>
  friend bool operator==(const ext_stack_allocator<T1, N1>& x,
                         const ext_stack_allocator<T1, N1>& y) noexcept;

 private:
  bool pointer_on_stack(T* ptr) const {
    return stack_ <= ptr && ptr < stack_ + N;
  }
};

template <class T, std::size_t N>
inline bool operator==(const ext_stack_allocator<T, N>& x,
                       const ext_stack_allocator<T, N>& y) noexcept {
  return x.stack_ == y.stack_ && x.free_ == y.free_;
}

template <class T, std::size_t N>
inline bool operator!=(const ext_stack_allocator<T, N>& x,
                       const ext_stack_allocator<T, N>& y) noexcept {
  return !(x == y);
}
}  // namespace detail
}  // namespace libint2

#endif  // INCLUDE_LIBINT2_UTIL_ARRAY_ADAPTOR_H_
