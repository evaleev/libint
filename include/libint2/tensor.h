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

#ifndef _libint2_include_tensor_h_
#define _libint2_include_tensor_h_

#include <libint2/util/cxxstd.h>
#if LIBINT2_CPLUSPLUS_STD < 2011
#error "The simple Libint API requires C++11 support"
#endif

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <numeric>
#include <vector>

namespace libint2 {
template <typename T>
struct Tensor {
 public:
  Tensor() = default;
  Tensor(const Tensor&) = default;
  Tensor(Tensor&&) = default;
  ~Tensor() = default;

  template <class... Dims>
  Tensor(Dims... dims) : dims_{std::forward<Dims>(dims)...} {
    strides_.resize(sizeof...(dims));
    // used in transform to compute strides
    struct stride {
      size_t value;
      stride() : value(1) {}
      size_t operator()(size_t dim) {
        auto old_value = value;
        value *= dim;
        return old_value;
      }
    };
    // row-major order of dimensions
    std::transform(dims_.rbegin(), dims_.rend(), strides_.rbegin(), stride());
    size_t size = strides_.size() == 0 ? 0 : strides_[0] * dims_[0];
    data_.resize(size);
  }

  T* data() { return &data_[0]; }
  const T* data() const { return static_cast<const T*>(this->data()); }
  template <class... MultiIndex>
  T* data(MultiIndex... index) {
    assert(sizeof...(MultiIndex) <= dims_.size());
    auto index_list = {std::forward<MultiIndex>(index)...};
    size_t ordinal = std::inner_product(index_list.begin(), index_list.end(),
                                        strides_.begin(), 0);
    return &data_[ordinal];
  }
  template <class... MultiIndex>
  const T* data(MultiIndex... index) const {
    return static_cast<const T*>(this->data());
  }

  template <class... MultiIndex>
  const T& operator()(MultiIndex... index) {
    assert(sizeof...(MultiIndex) == dims_.size());
    auto index_list = {std::forward<MultiIndex>(index)...};
    size_t ordinal = std::inner_product(index_list.begin(), index_list.end(),
                                        strides_.begin(), 0);
    return data_[ordinal];
  }

 private:
  std::vector<size_t> dims_;
  std::vector<size_t> strides_;
  std::vector<T> data_;

};  // class definition
}  // namespace libint2

#endif /* _libint2_include_tensor_h_*/
