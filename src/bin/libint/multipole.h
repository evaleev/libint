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

#ifndef _libint2_src_bin_libint_multipole_h_
#define _libint2_src_bin_libint_multipole_h_

#include <global_macros.h>
#include <hashable.h>
#include <polyconstr.h>
#include <smart_ptr.h>
#include <util_types.h>

#include <array>
#include <cassert>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

namespace libint2 {

/** Represents quantum numbers of cartesian multipole operator
 */
template <unsigned NDIM = 3>
class CartesianMultipoleQuanta
    : public Hashable<LIBINT2_UINT_LEAST64, ReferToKey> {
  static_assert(
      NDIM == 1 || NDIM == 3,
      "CartesianMultipoleQuanta<NDIM>: only NDIM=1 or NDIM=3 are supported");

 public:
  CartesianMultipoleQuanta() : valid_(true) {
    for (auto d = 0u; d != NDIM; ++d) n_[d] = 0u;
  }
  CartesianMultipoleQuanta(const CartesianMultipoleQuanta& other)
      : valid_(true) {
    std::copy(other.n_, other.n_ + NDIM, n_);
  }
  CartesianMultipoleQuanta& operator=(const CartesianMultipoleQuanta& other) {
    valid_ = other.valid_;
    std::copy(other.n_, other.n_ + NDIM, n_);
    return *this;
  }
  CartesianMultipoleQuanta& operator+=(const CartesianMultipoleQuanta& other) {
    assert(valid_);
    for (auto d = 0u; d != NDIM; ++d) n_[d] += other.n_[d];
    return *this;
  }
  CartesianMultipoleQuanta& operator-=(const CartesianMultipoleQuanta& other) {
    assert(valid_);
    for (auto d = 0u; d != NDIM; ++d) n_[d] -= other.n_[d];
    return *this;
  }
  ~CartesianMultipoleQuanta() = default;

  /// returns the number of quanta along xyz
  unsigned int operator[](unsigned int xyz) const {
    assert(xyz < NDIM);
    return n_[xyz];
  }
  /// Add c quanta along xyz.
  void inc(unsigned int xyz, unsigned int c = 1u) {
    assert(xyz < NDIM);
    assert(valid_);
    n_[xyz] += c;
  }
  /// Subtract c quanta along xyz. If impossible, invalidate the object, but do
  /// not change its quanta!
  void dec(unsigned int xyz, unsigned int c = 1u) {
    assert(xyz < NDIM);
    // assert(valid_);
    if (n_[xyz] >= c)
      n_[xyz] -= c;
    else
      valid_ = false;
  }
  /// Returns the sum of quantum numbers
  unsigned int norm() const { return std::accumulate(n_, n_ + NDIM, 0u); }
  /// norm() == 0
  bool zero() const { return norm() == 0; }
  /// Return false if this object is invalid
  bool valid() const { return valid_; }
  /// Implements Hashable<unsigned>::key()
  LIBINT2_UINT_LEAST64 key() const {
    if (NDIM == 3u) {
      unsigned nxy = n_[1] + n_[2];
      unsigned l = nxy + n_[0];
      LIBINT2_UINT_LEAST64 key = nxy * (nxy + 1) / 2 + n_[2];
      const auto result = key + key_l_offset.at(l);
      assert(result < max_key);
      return result;
    }
    if (NDIM == 1u) {
      const auto result = n_[0];
      assert(result < max_key);
      return result;
    }
    assert(false);
  }
  /// Return a compact label
  std::string label() const {
    char result[NDIM + 1];
    for (auto xyz = 0u; xyz < NDIM; ++xyz) result[xyz] = '0' + n_[xyz];
    result[NDIM] = '\0';
    return std::string(result);
  }

  /* ---------------
   * key computation
   * --------------- */
  constexpr static unsigned max_qn = LIBINT_CARTGAUSS_MAX_AM;

  // nkeys[L] is the number of possible CartesianMultipoleQuanta with \c norm()
  // == L \note for NDIM=1 nkeys[L] = 1 \note for NDIM=3 nkeys[L] = (L+1)(L+2)/2

  /// The range of keys is [0,max_key).
  /// \note for NDIM=3 the formula is easily derived by summing (L+1)(L+2)/2 up
  /// to \c max_qn
  const static unsigned max_key =
      NDIM == 3 ? (1 + max_qn) * (2 + max_qn) * (3 + max_qn) / 6 : (1 + max_qn);

  /// Print out the content
  void print(std::ostream& os = std::cout) const;

 protected:
  /// make this object invalid
  void invalidate() { valid_ = false; }

 private:
  unsigned int n_[NDIM];
  bool valid_;  // indicates valid/invalid state
  /// key_l_offset[L] is the number of all possible quanta combinations of order
  /// up to L \note key_l_offset[L] = sum k=[0,L) nkeys[L]
  static std::array<LIBINT2_UINT_LEAST64, CartesianMultipoleQuanta::max_qn + 1>
      key_l_offset;
};

// explicit instantiation declaration, definitions are multipole.cc
template <>
std::array<LIBINT2_UINT_LEAST64, CartesianMultipoleQuanta<1u>::max_qn + 1>
    CartesianMultipoleQuanta<1u>::key_l_offset;
template <>
std::array<LIBINT2_UINT_LEAST64, CartesianMultipoleQuanta<3u>::max_qn + 1>
    CartesianMultipoleQuanta<3u>::key_l_offset;

template <unsigned NDIM>
CartesianMultipoleQuanta<NDIM> operator-(
    const CartesianMultipoleQuanta<NDIM>& A,
    const CartesianMultipoleQuanta<NDIM>& B) {
  CartesianMultipoleQuanta<NDIM> Diff(A);
  for (unsigned int xyz = 0; xyz < 3; ++xyz) Diff.dec(xyz, B[xyz]);
  return Diff;
}

template <unsigned NDIM>
bool operator==(const CartesianMultipoleQuanta<NDIM>& A,
                const CartesianMultipoleQuanta<NDIM>& B) {
  for (unsigned d = 0; d != NDIM; ++d)
    if (A[d] != B[d]) return false;
  return true;
}

/// Return true if A is valid
template <unsigned NDIM>
inline bool exists(const CartesianMultipoleQuanta<NDIM>& A) {
  return A.valid();
}

////////////////////////////////////////////////////////////////////////////////////

/** Represents quantum numbers of \em real spherical multipole operator
 * defined in Eqs. 5 and 6 of J.M. Pérez-Jordá and W. Yang, J Chem Phys 104,
 * 8003 (1996). \f$ m \geq 0 \f$ corresponds to moments \f$ \mathcal{N}^+ \f$ ,
 * \f$ m < 0 \f$ corresponds to \f$ \mathcal{N}^- \f$ . To obtain the real solid
 * harmonics $C^m_l$ and $S^m_l$ defined in
 * https://en.wikipedia.org/wiki/Solid_harmonics multiply these harmonics by \f$
 * (-1)^m \sqrt{(2 - \delta_{m,0}) (l + |m|)! (l - |m|)!} \f$ .
 */
class SphericalMultipoleQuanta
    : public Hashable<LIBINT2_UINT_LEAST64, ReferToKey> {
 public:
  enum Sign { plus, minus };
  constexpr static unsigned max_qn = LIBINT_CARTGAUSS_MAX_AM;

  /// constructs an object in default (unusable) state
  SphericalMultipoleQuanta() : SphericalMultipoleQuanta(0, 1) {}
  /// constructs \f$ \mathcal{N}^{+}_{l,m} \f$ if \f$ m \geq 0 \f$, otherwise
  /// constructs \f$ \mathcal{N}^{-}_{l,m} \f$
  SphericalMultipoleQuanta(int l, int m)
      : SphericalMultipoleQuanta(l, m, m < 0 ? Sign::minus : Sign::plus) {}
  /// constructs \f$ \mathcal{N}^{\pm}_{l,m} \f$
  SphericalMultipoleQuanta(int l, int m, Sign sign)
      : l_(l),
        m_(std::abs(m)),
        sign_(sign),
        valid_(true),
        phase_(sign == Sign::plus && m < 0 ? -1 : 1) {
    if (l < 0) valid_ = false;
    if (m_ > l_) valid_ = false;
    // N^-_{0,0} = 0
    if (sign_ == Sign::minus && m_ == 0) valid_ = false;
  }

  int l() const {
    assert(valid_);
    return static_cast<int>(l_);
  }
  int m() const {
    assert(valid_);
    return static_cast<int>(m_);
  }
  Sign sign() const {
    assert(valid_);
    return sign_;
  }
  bool valid() const { return valid_; }
  int phase() const {
    assert(valid_);
    return phase_;
  }
  /// \f$ \mathcal{N}^{+}_{0,0} = 1 \f$
  bool is_precomputed() const {
    assert(valid_);
    return sign_ == Sign::plus && l_ == 0 && m_ == 0;
  }
  int value() const {
    assert(is_precomputed());
    return 1;
  }

  const static unsigned max_key = (1 + max_qn) * (1 + max_qn);

  /// Implements Hashable<unsigned>::key()
  LIBINT2_UINT_LEAST64 key() const {
    assert(valid_);
    const auto result = l_ * l_ + (sign_ == Sign::plus ? (l_ + m_) : (l_ - m_));
    assert(result < max_key);
    return result;
  }

 private:
  // N^-_{l,m} is symmetric is stored as N^-_{l,|m|}
  int l_;
  int m_;
  Sign sign_;
  bool valid_;
  int phase_;
};

/// Return true if A is valid
inline bool exists(const SphericalMultipoleQuanta& A) { return A.valid(); }

/// plus <-> minus
inline SphericalMultipoleQuanta::Sign flip(
    const SphericalMultipoleQuanta::Sign& s) {
  return s == SphericalMultipoleQuanta::Sign::minus
             ? SphericalMultipoleQuanta::Sign::plus
             : SphericalMultipoleQuanta::Sign::minus;
}

}  // namespace libint2

#endif
