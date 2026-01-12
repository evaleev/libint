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

#ifndef _libint2_include_cartesian_h_
#define _libint2_include_cartesian_h_

#include <libint2/cgshell_ordering.h>
#include <libint2/shell.h>

#include <array>

namespace libint2 {

namespace detail {

template <typename Real, std::size_t N>
struct scale; /* {
  void operator()(Real* intset, const std::array<std::pair<Real*,size_t>, N>&
coeffs);
}; */

template <typename Real>
struct scale<Real, 2> {
  inline void operator()(
      Real* intset, const std::array<std::pair<Real*, size_t>, 2>& coeffs) {
    auto* data = intset;
    for (auto f0 = 0ul; f0 != coeffs[0].second; ++f0) {
      for (auto f1 = 0ul; f1 != coeffs[1].second; ++f1) {
        const auto scalar01 = coeffs[0].first[f0] * coeffs[1].first[f1];
        *data *= scalar01;
        ++data;
      }
    }
  }
};

template <typename Real>
struct scale<Real, 4> {
  inline void operator()(
      Real* intset, const std::array<std::pair<Real*, size_t>, 4>& coeffs) {
    auto* data = intset;
    for (auto f0 = 0ul; f0 != coeffs[0].second; ++f0) {
      for (auto f1 = 0ul; f1 != coeffs[1].second; ++f1) {
        const auto scalar01 = coeffs[0].first[f0] * coeffs[1].first[f1];
        for (auto f2 = 0ul; f2 != coeffs[2].second; ++f2) {
          const auto scalar012 = scalar01 * coeffs[2].first[f2];
          for (auto f3 = 0ul; f3 != coeffs[3].second; ++f3) {
            const auto scalar0123 = scalar012 * coeffs[3].first[f3];
            *data *= scalar0123;
            ++data;
          }
        }
      }
    }
  }
};

// df_of_i_minux_1[i] = (i-1)!! , df_of_i_minux_1[0] = 1, df_of_i_minux_1[1] =
// 1, df_of_i_minux_1[2] = 1, df_of_i_minux_1[3] = 2, df_of_i_minux_1[4] = 3,
// etc.
template <typename Real>
inline std::vector<Real> make_df_of_i_minux_1(int imax) {
  std::vector<Real> df_of_i_minux_1(std::max(2, imax + 1));
  df_of_i_minux_1[0] = 1;
  df_of_i_minux_1[1] = 1;
  for (int i = 2; i <= imax; ++i) {
    df_of_i_minux_1[i] =
        (i - 1) * df_of_i_minux_1[i - 2];  // (i-1)!! = (i-1) * (i-3)!!
  }
  return df_of_i_minux_1;
}

template <typename Real>
inline std::vector<std::vector<Real>> make_cart_coeffs(int lmax) {
  static std::vector<Real> dfm1 =
      make_df_of_i_minux_1<Real>(2 * lmax);  // dfm1[i] = (i-1)!! , dfm1[0] = 1

  std::vector<std::vector<Real>> result(lmax + 1);
  for (int l = 0ul; l != lmax; ++l) {
    const auto cart_shell_size = (l + 1) * (l + 2) / 2;
    result[l].resize(cart_shell_size);
    int ixyz = 0;
    int ix, iy, iz;
    FOR_CART(ix, iy, iz, l)
    using std::sqrt;
    result[l][ixyz] = sqrt(
        dfm1.at(2 * l) / (dfm1.at(2 * ix) * dfm1.at(2 * iy) * dfm1.at(2 * iz)));
    ++ixyz;
    END_FOR_CART
  }

  return result;
}

}  // namespace detail

/// rescales cartesian Gaussians to convert from standard to uniform-normalized
/// convention
template <typename Real, std::size_t N>
inline void uniform_normalize_cartesian_shells(
    Real* intset, std::array<std::reference_wrapper<const Shell>, N> shells) {
  static std::vector<std::vector<Real>> cart_coeffs =
      detail::make_cart_coeffs<Real>(LIBINT_CARTGAUSS_MAX_AM);
  const auto max_shellsize_pure = 2 * LIBINT_CARTGAUSS_MAX_AM + 1;
  static std::vector<Real> pure_coeffs(max_shellsize_pure, Real(1));

  std::array<std::pair<Real*, size_t>, N> coeffs;
  for (auto c = 0u; c != N; ++c) {
    coeffs[c] =
        std::make_pair(shells[c].get().contr[0].pure
                           ? &pure_coeffs[0]
                           : &cart_coeffs[shells[c].get().contr[0].l][0],
                       shells[c].get().size());
  }

  detail::scale<Real, N>{}(intset, coeffs);
};

}  // namespace libint2

#endif  //_libint2_include_cartesian_h_
