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

#ifndef _libint2_src_bin_libint_bfset_h_
#define _libint2_src_bin_libint_bfset_h_

#include <contractable.h>
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

/** Set of basis functions. Sets must be constructable using
    std::shared_ptr<BFSet> or std::shared_ptr<ConstructablePolymorphically>.
*/
class BFSet : public ConstructablePolymorphically {
 public:
  virtual ~BFSet() {}
  virtual unsigned int num_bf() const = 0;
  virtual std::string label() const = 0;

 protected:
  BFSet() {}
};

/** Set of basis functions with incrementable/decrementable quantum numbers.
    Sets must be constructable using std::shared_ptr<BFSet> or
   std::shared_ptr<ConstructablePolymorphically>.

    Call to dec() may invalidate the object. No further
    modification of such object's state is possible.
    Assignment of such object to another preserves the invalid state,
    but copy construction of a valid object is possible.
*/
class IncableBFSet : public BFSet {
 public:
  virtual ~IncableBFSet() {}

  /// Add c quanta along xyz.
  virtual void inc(unsigned int xyz, unsigned int c = 1u) = 0;
  /// Subtract c quanta along xyz. If impossible, invalidate the object, but do
  /// not change its quanta!
  virtual void dec(unsigned int xyz, unsigned int c = 1u) = 0;
  /// Returns the norm of the quantum numbers
  virtual unsigned int norm() const = 0;
  /// norm() == 0
  bool zero() const { return norm() == 0; }
  /// Return false if this object is invalid
  bool valid() const { return valid_; }

 protected:
  IncableBFSet() : valid_(true) {}

  /// make this object invalid
  void invalidate() { valid_ = false; }

 private:
  bool valid_;
};

/// BF with unit quantum along X. F must behave like IncableBFSet
template <typename F>
F unit(unsigned int X) {
  F tmp;
  tmp.inc(X, 1);
  return tmp;
}
/// Return true if A is valid
inline bool exists(const IncableBFSet& A) { return A.valid(); }

/** Represents cartesian derivatives of atom-centered basis functions
 */
template <unsigned NDIM = 3>
class OriginDerivative : public Hashable<LIBINT2_UINT_LEAST64, ReferToKey> {
  static_assert(NDIM == 1 || NDIM == 3,
                "OriginDerivative<NDIM>: only NDIM=1 or NDIM=3 are supported");

 public:
  OriginDerivative() : valid_(true) {
    for (auto d = 0u; d != NDIM; ++d) d_[d] = 0u;
  }
  OriginDerivative(const OriginDerivative& other) : valid_(true) {
    std::copy(other.d_, other.d_ + NDIM, d_);
  }
  OriginDerivative& operator=(const OriginDerivative& other) {
    valid_ = other.valid_;
    std::copy(other.d_, other.d_ + NDIM, d_);
    return *this;
  }
  OriginDerivative& operator+=(const OriginDerivative& other) {
    assert(valid_);
    for (auto d = 0u; d != NDIM; ++d) d_[d] += other.d_[d];
    return *this;
  }
  OriginDerivative& operator-=(const OriginDerivative& other) {
    assert(valid_);
    for (auto d = 0u; d != NDIM; ++d) d_[d] -= other.d_[d];
    return *this;
  }
  ~OriginDerivative() {
    static_assert(NDIM == 3u || NDIM == 1u,
                  "OriginDerivative with NDIM=1,3 are implemented");
  }

  /// returns the number of quanta along xyz
  unsigned int d(unsigned int xyz) const {
    assert(xyz < NDIM);
    return d_[xyz];
  }
  /// returns the number of quanta along xyz
  unsigned int operator[](unsigned int xyz) const { return this->d(xyz); }
  /// Add c quanta along xyz.
  void inc(unsigned int xyz, unsigned int c = 1u) {
    assert(xyz < NDIM);
    assert(valid_);
    d_[xyz] += c;
  }
  /// Subtract c quanta along xyz. If impossible, invalidate the object, but do
  /// not change its quanta!
  void dec(unsigned int xyz, unsigned int c = 1u) {
    assert(xyz < NDIM);
    // assert(valid_);
    if (d_[xyz] >= c)
      d_[xyz] -= c;
    else
      valid_ = false;
  }
  /// Returns the norm of the quantum numbers
  unsigned int norm() const { return std::accumulate(d_, d_ + NDIM, 0u); }
  /// norm() == 0
  bool zero() const { return norm() == 0; }
  /// Return false if this object is invalid
  bool valid() const { return valid_; }
  /// Implements Hashable<unsigned>::key()
  LIBINT2_UINT_LEAST64 key() const override {
    if (NDIM == 3u) {
      unsigned nxy = d_[1] + d_[2];
      unsigned l = nxy + d_[0];
      LIBINT2_UINT_LEAST64 key = nxy * (nxy + 1) / 2 + d_[2];
      const auto result = key + key_l_offset.at(l);
      assert(result < max_key);
      return result;
    }
    if (NDIM == 1u) {
      const auto result = d_[0];
      assert(result < max_key);
      return result;
    }
    assert(false);
  }
  /// Return a compact label
  std::string label() const {
    char result[NDIM + 1];
    for (auto xyz = 0u; xyz < NDIM; ++xyz) result[xyz] = '0' + d_[xyz];
    result[NDIM] = '\0';
    return std::string(result);
  }

  /* ---------------
   * key computation
   * --------------- */
  const static unsigned max_deriv = 4;

  // nkeys_l[L] is the number of possible OriginDerivative's with \c norm() == L
  // \note for NDIM=1 nkeys_l[L] = 1
  // \note for NDIM=2 nkeys_l[L] = (L+1)
  // \note for NDIM=3 nkeys_l[L] = (L+1)(L+2)/2
  // static std::array<LIBINT2_UINT_LEAST64, OriginDerivative::max_deriv+1>
  // nkeys;

  /// The range of keys is [0,max_key).
  /// \note for NDIM=3 the formula is easily derived by summing (L+1)(L+2)/2 up
  /// to \c max_deriv
  const static unsigned max_key =
      NDIM == 3 ? (1 + max_deriv) * (2 + max_deriv) * (3 + max_deriv) / 6
                : (1 + max_deriv);

  /// Print out the content
  void print(std::ostream& os = std::cout) const;

 protected:
  /// make this object invalid
  void invalidate() { valid_ = false; }

 private:
  unsigned int d_[NDIM];
  bool valid_;  // indicates valid/invalid state
  /// key_l_offset[L] is the number of all possible derivatives of order up to L
  /// \note key_l_offset[L] = sum k=[0,L) nkeys[L]
  static std::array<LIBINT2_UINT_LEAST64, OriginDerivative::max_deriv + 1>
      key_l_offset;
};

// explicit instantiation declaration, definitions are gauss.cc
template <>
std::array<LIBINT2_UINT_LEAST64, OriginDerivative<1u>::max_deriv + 1>
    OriginDerivative<1u>::key_l_offset;
template <>
std::array<LIBINT2_UINT_LEAST64, OriginDerivative<3u>::max_deriv + 1>
    OriginDerivative<3u>::key_l_offset;

template <unsigned NDIM>
OriginDerivative<NDIM> operator-(const OriginDerivative<NDIM>& A,
                                 const OriginDerivative<NDIM>& B) {
  OriginDerivative<NDIM> Diff(A);
  for (unsigned int xyz = 0; xyz < 3; ++xyz) Diff.dec(xyz, B.d(xyz));
  return Diff;
}

template <unsigned NDIM>
bool operator==(const OriginDerivative<NDIM>& A,
                const OriginDerivative<NDIM>& B) {
  for (unsigned d = 0; d != NDIM; ++d)
    if (A.d(d) != B.d(d)) return false;
  return true;
}

/// Return true if A is valid
template <unsigned NDIM>
inline bool exists(const OriginDerivative<NDIM>& A) {
  return A.valid();
}

class CGF;  // forward declaration of CGF

/// 3D Cartesian Gaussian Shell
class CGShell : public IncableBFSet,
                public Hashable<LIBINT2_UINT_LEAST64, ReferToKey>,
                public Contractable<CGShell> {
  unsigned int qn_[1];
  OriginDerivative<3> deriv_;
  bool pure_sh_;  //< if true, assumed to contain solid harmonics with quantum
                  // number qn_[0] only
  /** if true, this is a unit shell (zero-exponent Gaussian) */
  bool unit_;

  friend CGShell operator+(const CGShell& A, const CGShell& B);
  friend CGShell operator-(const CGShell& A, const CGShell& B);

 public:
  /// As far as SetIterator is concerned, CGShell is a set of CGFs
  typedef CGF iter_type;
  typedef IncableBFSet parent_type;

  /// Default constructor creates an s-type shell
  CGShell();
  CGShell(unsigned int qn, bool pure_sh = false);
  CGShell(const CGShell&);
  virtual ~CGShell();
  CGShell& operator=(const CGShell&);

  const OriginDerivative<3u>& deriv() const { return deriv_; }
  OriginDerivative<3u>& deriv() { return deriv_; }

  /// Return a compact label
  std::string label() const override;
  /// Returns the number of basis functions in the set
  unsigned int num_bf() const override {
    return (qn_[0] + 1) * (qn_[0] + 2) / 2;
  }
  /// Returns the angular momentum
  unsigned int qn(unsigned int xyz = 0) const { return qn_[0]; }
  /// returns the angular momentum
  unsigned int operator[](unsigned int xyz) const { return this->qn(xyz); }

  /// Comparison operator
  bool operator==(const CGShell&) const;

  /// contains only solid harmonics with the same quantum number as this shell?
  /// (this may permit simplified RR to be used -- obviously must transform to
  /// solid harmonics later)
  bool pure_sh() const { return pure_sh_; }
  /// @param p if true, will assume to contain only solid harmonics of the same
  /// quantum number as this shell
  void pure_sh(bool p) { pure_sh_ = p; }

  /// Implementation of IncableBFSet::inc().
  void inc(unsigned int xyz, unsigned int c = 1u) override;
  /// Implementation of IncableBFSet::dec().
  void dec(unsigned int xyz, unsigned int c = 1u) override;
  /// Implements IncableBFSet::norm()
  unsigned int norm() const override;
  /// Implements Hashable<LIBINT2_UINT_LEAST64>::key()
  LIBINT2_UINT_LEAST64 key() const override {
    if (is_unit()) return max_key - 1;
    const LIBINT2_UINT_LEAST64 result =
        ((deriv().key() * 2 + (contracted() ? 1 : 0)) * (max_qn + 1) + qn_[0]) *
            2 +
        (pure_sh() ? 1 : 0);
    assert(result < max_key - 1);
    return result;
  }
  const static LIBINT2_UINT_LEAST64 max_qn = LIBINT_CARTGAUSS_MAX_AM;
  /** The range of keys is [0,max_key]
      deriv_key_range = OriginDerivative<3u>::max_key
      contracted = 2 (yes or no)
      qn_range = max_qn + 1
      puresh_key_range = 2
      +1 to account for the unit shell
    */
  const static LIBINT2_UINT_LEAST64 max_key =
      OriginDerivative<3u>::max_key * 2 * (max_qn + 1) * 2 + 1;

  /// Print out the content
  void print(std::ostream& os = std::cout) const;

  /// returns the unit shell (exponent=0, am=0, indicated by pure_sh=true)
  static CGShell unit();
  bool is_unit() const { return unit_; }
};

CGShell operator+(const CGShell& A, const CGShell& B);
CGShell operator-(const CGShell& A, const CGShell& B);

/// 3D Cartesian Gaussian Function
class CGF : public IncableBFSet,
            public Hashable<LIBINT2_UINT_LEAST64, ComputeKey>,
            public Contractable<CGF> {
  unsigned int qn_[3];
  OriginDerivative<3u> deriv_;
  bool pure_sh_;  //< if true, assumed to contain solid harmonics with quantum
                  // number qn_[0] only
  bool unit_;     //< if true, this is a unit Gaussian (exponent = 0)

  friend CGF operator+(const CGF& A, const CGF& B);
  friend CGF operator-(const CGF& A, const CGF& B);

 public:
  /// As far as SetIterator is concerned, CGF is a set of one CGF
  typedef CGF iter_type;
  typedef IncableBFSet parent_type;
  /// How to return key
  // typedef typename Hashable<unsigned,ComputeKey>::KeyReturnType
  // KeyReturnType;

  /// Default constructor makes an s-type Gaussian
  CGF();
  CGF(unsigned int qn[3], bool pure_sh = false);
  CGF(const CGF&);
  explicit CGF(const ConstructablePolymorphically&);
  virtual ~CGF();
  /// assignment
  CGF& operator=(const CGF&);

  const OriginDerivative<3u>& deriv() const { return deriv_; }
  OriginDerivative<3u>& deriv() { return deriv_; }

  /// Return a compact label
  std::string label() const override;
  /// Returns the number of basis functions in the set (always 1)
  unsigned int num_bf() const override { return 1; }
  /// Returns the quantum number along \c axis
  unsigned int qn(unsigned int axis) const;
  unsigned int operator[](unsigned int axis) const { return qn(axis); }

  /// contains only solid harmonics with the same quantum number as this shell?
  /// (this may permit simplified RR to be used -- obviously must transform to
  /// solid harmonics later)
  bool pure_sh() const { return pure_sh_; }
  /// @param p if true, will assume to contain only solid harmonics of the same
  /// quantum number as this shell
  void pure_sh(bool p) { pure_sh_ = p; }

  /// Comparison operator
  bool operator==(const CGF&) const;

  /// Implementation of IncableBFSet::inc().
  void inc(unsigned int xyz, unsigned int c = 1u) override;
  /// Implementation of IncableBFSet::dec().
  void dec(unsigned int xyz, unsigned int c = 1u) override;
  /// Implements IncableBFSet::norm()
  unsigned int norm() const override;
  /// Implements Hashable<LIBINT2_UINT_LEAST64>::key()
  LIBINT2_UINT_LEAST64 key() const override {
    if (is_unit()) return max_key - 1;
    unsigned nxy = qn_[1] + qn_[2];
    unsigned l = nxy + qn_[0];
    LIBINT2_UINT_LEAST64 key = nxy * (nxy + 1) / 2 + qn_[2];
    const LIBINT2_UINT_LEAST64 result =
        ((deriv().key() * 2 + (contracted() ? 1 : 0)) * max_num_qn + key +
         key_l_offset.at(l)) *
            2 +
        (pure_sh() ? 1 : 0);
    if (result >= max_key - 1) {
      this->print(std::cout);
      std::cout << "result,max_key-1 = " << result << "," << max_key - 1
                << std::endl;
      assert(result < max_key - 1);
    }
    return result;
  }
  /// The range of keys is [0,max_key). The formula is easily derived by summing
  /// (L+1)(L+2)/2 up to CGShell::max_key The factor of 2 to account for
  /// contracted vs. uncontracted basis functions The factor of
  /// OriginDerivative::max_key to account for derivatives
  const static LIBINT2_UINT_LEAST64 max_num_qn =
      ((1 + (CGShell::max_qn + 1)) * (2 + (CGShell::max_qn + 1)) *
       (3 + (CGShell::max_qn + 1)) / 6);
  // deriv_key_range = OriginDerivative<3u>::max_key
  // contracted = 2 (yes or no)
  // qn_range = max_num_qn
  // puresh_key_range = 2
  // +1 to account for unit function
  const static LIBINT2_UINT_LEAST64 max_key =
      OriginDerivative<3u>::max_key * 2ul * max_num_qn * 2ul + 1;

  /// Print out the content
  void print(std::ostream& os = std::cout) const;

  /// returns the unit shell (exponent=0, am=0, indicated by unit_=true)
  static CGF unit();
  bool is_unit() const { return unit_; }

 private:
  /// key_l_offset[L] is the number of all possible CGFs with angular momentum
  /// less than L
  static std::array<LIBINT2_UINT_LEAST64, CGShell::max_qn + 1> key_l_offset;
};

CGF operator+(const CGF& A, const CGF& B);
CGF operator-(const CGF& A, const CGF& B);

#if 1
/// Cartesian components of 3D CGF = 1D CGF
/// @note reference to particular cartesian axis embedded in type, so that for
/// example axis-dependent RRs can correctly infer geometric constants
template <CartesianAxis Axis>
class CGF1d : public IncableBFSet,
              public Hashable<LIBINT2_UINT_LEAST64, ComputeKey>,
              public Contractable<CGF1d<Axis> > {
  unsigned int qn_[1];
  OriginDerivative<1u> deriv_;
  bool unit_;  //< if true, this is a unit Gaussian (exponent = 0)

 public:
  static constexpr auto axis = Axis;

  /// As far as SetIterator is concerned, CGF1d is a set of one CGF1d
  typedef CGF1d iter_type;
  typedef IncableBFSet parent_type;

  /// Default constructor makes an qn=0 Gaussian
  CGF1d() : unit_(false) { qn_[0] = 0; }
  CGF1d(unsigned int qn) : unit_(false) { qn_[0] = qn; }
  CGF1d(unsigned int qn[1]) : unit_(false) { qn_[0] = qn[0]; }
  CGF1d(const CGF1d& source)
      : Contractable<CGF1d>(source),
        deriv_(source.deriv_),
        unit_(source.unit_) {
    qn_[0] = source.qn_[0];
  }
  explicit CGF1d(const ConstructablePolymorphically& sptr)
      : Contractable<CGF1d>(dynamic_cast<const CGF1d&>(sptr)) {
    const CGF1d& sptr_cast = dynamic_cast<const CGF1d&>(sptr);
    qn_[0] = sptr_cast.qn_[0];
    deriv_ = sptr_cast.deriv_;
    unit_ = sptr_cast.unit_;
  }
  virtual ~CGF1d() {}

  /// assignment
  CGF1d& operator=(const CGF1d& source) {
    qn_[0] = source.qn_[0];
    deriv_ = source.deriv_;
    unit_ = source.unit_;
    Contractable<CGF1d>::operator=(source);
    if (!source.valid()) invalidate();
    return *this;
  }

  CGF1d operator+(const CGF1d& B) const {
    // assert(this->is_unit() == false && B.is_unit() == false);
    CGF1d<Axis> Sum(*this);
    Sum.inc(0, B.qn(0));
    Sum.deriv_ += B.deriv_;
    return Sum;
  }
  CGF1d operator-(const CGF1d& B) const {
    // assert(A.is_unit() == false && B.is_unit() == false);
    CGF1d Diff(*this);
    Diff.dec(0, B.qn(0));
    Diff.deriv_ -= B.deriv_;
    return Diff;
  }

  const OriginDerivative<1u>& deriv() const { return deriv_; }
  OriginDerivative<1u>& deriv() { return deriv_; }

  /// Return a compact label
  std::string label() const override {
    // unit *functions* are treated as regular qn-0 functions so that
    // (00|00)^(m) = (unit 0|00)^(m)
    std::ostringstream oss;
    oss << to_string(Axis) << qn_[0];
    if (!deriv_.zero()) oss << "_" << deriv_.label();

    // I don't handle labels of contracted CGF1d because I don't think I need
    // them make sure just in case
    assert(!this->contracted());

    return oss.str();
  }

  /// Returns the number of basis functions in the set (always 1)
  unsigned int num_bf() const override { return 1; }
  /// Returns the quantum number (what used to be "angular momentum")
  unsigned int qn(unsigned int dir = 0) const {
    assert(dir == 0);
    return qn_[0];
  }
  unsigned int operator[](unsigned int dir) const { return this->qn(dir); }

  /// Comparison operator
  bool operator==(const CGF1d& a) const {
    return (qn_[0] == a.qn_[0] && this->contracted() == a.contracted() &&
            deriv_ == a.deriv_ && unit_ == a.unit_);
  }

  /// Implementation of IncableBFSet::inc().
  void inc(unsigned int dir, unsigned int c = 1u) override {
    assert(!is_unit());
    assert(dir == 0);
    if (valid()) qn_[0] += c;
  }
  /// Implementation of IncableBFSet::dec().
  void dec(unsigned int dir, unsigned int c = 1u) override {
    if (is_unit()) {
      invalidate();
      return;
    }
    assert(dir == 0);
    if (valid()) {
      if (qn_[0] < c) {
        invalidate();
        return;
      }
      qn_[0] -= c;
    }
  }
  /// Implements IncableBFSet::norm()
  unsigned int norm() const override { return qn_[0]; }
  /// Implements Hashable<LIBINT2_UINT_LEAST64>::key()
  LIBINT2_UINT_LEAST64 key() const override {
    if (is_unit()) return max_key - 1;
    const LIBINT2_UINT_LEAST64 result =
        (deriv().key() * 2ul + (this->contracted() ? 1ul : 0ul)) * max_num_qn +
        qn_[0];
    if (result >= max_key - 1) {
      this->print(std::cout);
      std::cout << "result,max_key-1 = " << result << "," << max_key - 1
                << std::endl;
      assert(result < max_key - 1);
    }
    return result;
  }
  /// The range of keys is [0,max_key). The formula is easily derived by summing
  /// (L+1)(L+2)/2 up to CGShell::max_key The factor of 2 to account for
  /// contracted vs. uncontracted basis functions The factor of
  /// OriginDerivative::max_key to account for derivatives
  const static LIBINT2_UINT_LEAST64 max_num_qn = CGShell::max_qn + 1;
  // deriv_key_range = OriginDerivative<1u>::max_key
  // contracted = 2 (yes or no)
  // qn_range = max_num_qn
  // +1 to account for unit function
  const static LIBINT2_UINT_LEAST64 max_key =
      OriginDerivative<1u>::max_key * OriginDerivative<1u>::max_key *
          max_num_qn +
      1;

  /// Print out the content
  void print(std::ostream& os = std::cout) const {
    os << "CGF1d<" << to_string(Axis) << ">: " << label() << std::endl;
  }

  /// returns the unit shell (exponent=0, am=0, indicated by unit_=true)
  static CGF1d unit() {
    CGF1d result;
    result.unit_ = true;
    result.uncontract();
    return result;
  }
  bool is_unit() const { return unit_; }

 private:
  /// key_l_offset[L] is the number of all possible CGF1d's with quantum number
  /// less than L
  static std::array<LIBINT2_UINT_LEAST64, CGShell::max_qn + 1> key_l_offset;
};

//  template <CartesianAxis Axis>
//  inline CGF1d<Axis> operator+(const CGF1d<Axis>& A, const CGF1d<Axis>& B) {
//    assert(A.is_unit() == false && B.is_unit() == false);
//    CGF1d<Axis> Sum(A);
//    Sum.inc(0,B.qn(0));
//    Sum.deriv_ += B.deriv_;
//    return Sum;
//  }
//  template <CartesianAxis Axis>
//  inline CGF1d<Axis> operator-(const CGF1d<Axis>& A, const CGF1d<Axis>& B) {
//    //assert(A.is_unit() == false && B.is_unit() == false);
//    CGF1d<Axis> Diff(A);
//    Diff.dec(0,B.qn(0));
//    Diff.deriv_ -= B.deriv_;
//
//    return Diff;
//  }

/// a "shell" of 1D CGFs with quantum number L is a set of 1D CGFs with quantum
/// numbers 0 .. L
///
/// @note This is very different from a CGShell which consists of CGFs with same
/// "norm", distributed
///       differently between axes. The notion of 1d shell is still useful
///       because we want to compute integrals over all functions in the shell
///       at once.
/// @note Just like with CGF1d, the axis is embedded into the type
template <CartesianAxis Axis>
class CGShell1d : public IncableBFSet,
                  public Hashable<LIBINT2_UINT_LEAST64, ComputeKey>,
                  public Contractable<CGShell1d<Axis> > {
  unsigned int qn_[1];
  OriginDerivative<1u> deriv_;
  bool unit_;  //< if true, this is a unit Gaussian (exponent = 0)

 public:
  static constexpr auto axis = Axis;

  /// CGShell1d is a set CGF1d's
  typedef CGF1d<Axis> iter_type;
  typedef IncableBFSet parent_type;

  /// Default constructor makes a qn=0 shell
  CGShell1d() : unit_(false) { qn_[0] = 0; }
  CGShell1d(unsigned int qn) : unit_(false) { qn_[0] = qn; }
  CGShell1d(unsigned int qn[1]) : unit_(false) { qn_[0] = qn[0]; }
  CGShell1d(const CGShell1d& source)
      : Contractable<CGShell1d>(source),
        deriv_(source.deriv_),
        unit_(source.unit_) {
    qn_[0] = source.qn_[0];
  }
  virtual ~CGShell1d() {}

  /// assignment
  CGShell1d& operator=(const CGShell1d& source) {
    qn_[0] = source.qn_[0];
    deriv_ = source.deriv_;
    unit_ = source.unit_;
    Contractable<CGShell1d>::operator=(source);
    if (!source.valid()) invalidate();
    return *this;
  }

  const OriginDerivative<1u>& deriv() const { return deriv_; }
  OriginDerivative<1u>& deriv() { return deriv_; }

  /// Return a compact label
  std::string label() const override {
    // unit *functions* are treated as regular qn-0 functions so that
    // (00|00)^(m) = (unit 0|00)^(m)
    std::ostringstream oss;
    auto axis_label = to_string(Axis);
    axis_label[0] = std::toupper(axis_label[0]);
    oss << axis_label << qn_[0];
    if (!deriv_.zero()) oss << "_" << deriv_.label();

    // I don't handle labels of contracted CGF1d because I don't think I need
    // them make sure just in case
    assert(!this->contracted());

    return oss.str();
  }

  /// Returns the number of basis functions in the set (always 1)
  unsigned int num_bf() const override { return qn_[0] + 1; }
  /// Returns the quantum number (what used to be "angular momentum")
  unsigned int qn(unsigned int dir = 0) const {
    assert(dir == 0);
    return qn_[0];
  }

  /// Comparison operator
  bool operator==(const CGShell1d& a) const {
    return (qn_[0] == a.qn_[0] && this->contracted() == a.contracted() &&
            deriv_ == a.deriv_ && unit_ == a.unit_);
  }

  /// Implementation of IncableBFSet::inc().
  void inc(unsigned int dir, unsigned int c = 1u) override { assert(false); }
  /// Implementation of IncableBFSet::dec().
  void dec(unsigned int dir, unsigned int c = 1u) override { assert(false); }
  /// Implements IncableBFSet::norm()
  unsigned int norm() const override { return qn_[0]; }
  /// Implements Hashable<LIBINT2_UINT_LEAST64>::key()
  LIBINT2_UINT_LEAST64 key() const override {
    if (is_unit()) return max_key - 1;
    const LIBINT2_UINT_LEAST64 result =
        (deriv().key() * 2ul + (this->contracted() ? 1ul : 0ul)) * max_num_qn +
        qn_[0];
    if (result >= max_key - 1) {
      this->print(std::cout);
      std::cout << "result,max_key-1 = " << result << "," << max_key - 1
                << std::endl;
      assert(result < max_key - 1);
    }
    return result;
  }
  /// The range of keys is [0,max_key). The formula is easily derived by summing
  /// (L+1)(L+2)/2 up to CGShell::max_key The factor of 2 to account for
  /// contracted vs. uncontracted basis functions The factor of
  /// OriginDerivative::max_key to account for derivatives
  const static LIBINT2_UINT_LEAST64 max_num_qn = CGShell::max_qn + 1;
  // deriv_key_range = OriginDerivative<1u>::max_key
  // contracted = 2 (yes or no)
  // qn_range = max_num_qn
  // +1 to account for unit function
  const static LIBINT2_UINT_LEAST64 max_key =
      OriginDerivative<1u>::max_key * OriginDerivative<1u>::max_key *
          max_num_qn +
      1;

  /// Print out the content
  void print(std::ostream& os = std::cout) const {
    os << "CGShell1d<" << to_string(Axis) << ">: " << label() << std::endl;
  }

  /// returns the unit shell (exponent=0, am=0, indicated by unit_=true)
  static CGShell1d unit() {
    CGShell1d result;
    result.unit_ = true;
    result.uncontract();
    return result;
  }
  bool is_unit() const { return unit_; }

 private:
  /// key_l_offset[L] is the number of all possible CGF1d's with quantum number
  /// less than L
  static std::array<LIBINT2_UINT_LEAST64, CGShell::max_qn + 1> key_l_offset;
};

#endif

#if 0
  class SHGF; // forward declaration

  /// Solid-Harmonic Gaussian Shell
  class SHGShell : public IncableBFSet, public Hashable<unsigned,ReferToKey>,
                   public Contractable<SHGShell> {

    unsigned int qn_[1];
    OriginDerivative deriv_;

    friend SHGShell operator+(const SHGShell& A, const SHGShell& B);
    friend SHGShell operator-(const SHGShell& A, const SHGShell& B);

  public:
    /// As far as SetIterator is concerned, SHGShell is a set of SHGFs
    typedef SHGF iter_type;
    typedef IncableBFSet parent_type;

    /// Default constructor creates an s-type shell
    SHGShell();
    SHGShell(unsigned int qn);
    SHGShell(const SHGShell&);
    virtual ~SHGShell();
    SHGShell& operator=(const SHGShell&);

    const OriginDerivative& deriv() const { return deriv_; }
    OriginDerivative& deriv() { return deriv_; }

    /// Return a compact label
    std::string label() const;
    /// Returns the number of basis functions in the set
    unsigned int num_bf() const { return 2*qn_[0]+1; };
    /// Returns the angular momentum
    unsigned int qn(unsigned int m=0) const { return qn_[0]; }

    /// Comparison operator
    bool operator==(const SHGShell&) const;

    /// Implementation of IncableBFSet::inc().
    void inc(unsigned int xyz, unsigned int c = 1u);
    /// Implementation of IncableBFSet::dec().
    void dec(unsigned int xyz, unsigned int c = 1u);
    /// Implements IncableBFSet::norm()
    unsigned int norm() const;
    /// Implements Hashable<unsigned>::key()
    unsigned key() const { return (deriv().key() * 2 + (contracted() ? 1 : 0)) * (max_qn+1) + qn_[0]; }
    const static unsigned max_qn = LIBINT_CARTGAUSS_MAX_AM;
    // The range of keys is [0,max_key]
    const static unsigned max_key = 2 * (max_qn + 1) * OriginDerivative::max_key * 2; // deriv_key_range = 2
                                                                                      // qn_range = max_qn + 1
                                                                                      // puresh_key_range = 2

    /// Print out the content
    void print(std::ostream& os = std::cout) const;

  };

  SHGShell operator+(const SHGShell& A, const SHGShell& B);
  SHGShell operator-(const SHGShell& A, const SHGShell& B);

  /// Solid-Harmonic Gaussian Function
  class SHGF : public IncableBFSet, public Hashable<unsigned,ComputeKey>,
               public Contractable<SHGF> {

    unsigned int qn_[3];
    OriginDerivative deriv_;

    friend SHGF operator+(const SHGF& A, const SHGF& B);
    friend SHGF operator-(const SHGF& A, const SHGF& B);

  public:
    /// As far as SetIterator is concerned, SHGF is a set of one SHGF
    typedef SHGF iter_type;
    typedef IncableBFSet parent_type;
    /// How to return key
    //typedef typename Hashable<unsigned,ComputeKey>::KeyReturnType KeyReturnType;

    /// Default constructor makes an s-type Gaussian
    SHGF();
    SHGF(unsigned int qn[3]);
    SHGF(const SHGF&);
    SHGF(const ConstructablePolymorphically&);
    virtual ~SHGF();
    /// assignment
    SHGF& operator=(const SHGF&);

    const OriginDerivative& deriv() const { return deriv_; }
    OriginDerivative& deriv() { return deriv_; }

    /// Return a compact label
    std::string label() const;
    /// Returns the number of basis functions in the set (always 1)
    unsigned int num_bf() const { return 1; };
    /// Returns the angular momentum
    unsigned int qn(unsigned int xyz) const;

    /// Comparison operator
    bool operator==(const SHGF&) const;

    /// Implementation of IncableBFSet::inc().
    void inc(unsigned int xyz, unsigned int c = 1u);
    /// Implementation of IncableBFSet::dec().
    void dec(unsigned int xyz, unsigned int c = 1u);
    /// Implements IncableBFSet::norm()
    unsigned int norm() const;
    /// Implements Hashable<unsigned>::key()
    unsigned key() const {
      unsigned nxy = qn_[1] + qn_[2];
      unsigned l = nxy + qn_[0];
      unsigned key = nxy*(nxy+1)/2 + qn_[2];
      return ( deriv().key() * 2 + (contracted() ? 1 : 0)) * max_num_qn + key + key_l_offset.at(l);
    }
    /// The range of keys is [0,max_key). The formula is easily derived by summing (L+1)(L+2)/2 up to SHGShell::max_key
    /// The factor of 2 to account for contracted vs. uncontracted basis functions
    /// The factor of OriginDerivative::max_key to account for derivatives
    const static unsigned max_num_qn = ((1 + (SHGShell::max_qn+1)) * (2 + (SHGShell::max_qn+1)) * (3 + (SHGShell::max_qn+1)) /6);
    const static unsigned max_key = 2 * OriginDerivative::max_key * max_num_qn;

    /// Print out the content
    void print(std::ostream& os = std::cout) const;

  private:
    /// key_l_offset[L] is the number of all possible SGGFs with angular momentum less than L
    static unsigned key_l_offset[SHGShell::max_key+2];
  };

  SHGF operator+(const SHGF& A, const SHGF& B);
  SHGF operator-(const SHGF& A, const SHGF& B);
#endif

/**
   TrivialBFSet<T> defines static member result, which is true if T
   is a basis function set consisting of 1 function
*/
template <class T>
struct TrivialBFSet;
template <>
struct TrivialBFSet<CGShell> {
  static const bool result = false;
};
template <>
struct TrivialBFSet<CGF> {
  static const bool result = true;
};
template <CartesianAxis Axis>
struct TrivialBFSet<CGShell1d<Axis> > {
  static const bool result = false;
};
template <CartesianAxis Axis>
struct TrivialBFSet<CGF1d<Axis> > {
  static const bool result = true;
};
#if 0
  template <>
    struct TrivialBFSet<SHGShell> {
      static const bool result = false;
    };
  template <>
    struct TrivialBFSet<SHGF> {
      static const bool result = true;
    };
#endif

};  // namespace libint2

#endif
