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
#include <ctype.h>
#include <default_params.h>
#include <exception.h>
#include <multipole.h>
#include <string.h>

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace std;
using namespace libint2;

namespace libint2 {
extern LIBINT2_UINT_LEAST64 cgshell_key_l_offset_[];
}

namespace {
template <typename T, std::size_t N>
std::array<T, N> make_std_array(T* data) {
  std::array<T, N> result;
  std::copy(data, data + N, result.begin());
  return result;
}
template <typename T, std::size_t N>
std::array<T, N> make_std_array(const std::initializer_list<T>& data) {
  assert(data.size() >= N);
  std::array<T, N> result;
  std::copy(data.begin(), data.begin() + N, result.begin());
  return result;
}
};  // namespace

template <>
std::array<LIBINT2_UINT_LEAST64, CartesianMultipoleQuanta<1u>::max_qn + 1>
    CartesianMultipoleQuanta<1u>::key_l_offset(
        make_std_array<LIBINT2_UINT_LEAST64,
                       CartesianMultipoleQuanta::max_qn + 1>(
            {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
             12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
             24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34}));
template <>
std::array<LIBINT2_UINT_LEAST64, CartesianMultipoleQuanta<3u>::max_qn + 1>
    CartesianMultipoleQuanta<3u>::key_l_offset(
        make_std_array<LIBINT2_UINT_LEAST64, CartesianMultipoleQuanta::max_qn +
                                                 1>(cgshell_key_l_offset_));
