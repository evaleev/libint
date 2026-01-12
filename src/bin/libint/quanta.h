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

#ifndef _libint2_src_bin_libint_quanta_h_
#define _libint2_src_bin_libint_quanta_h_

#include <global_macros.h>
#include <iter.h>
#include <smart_ptr.h>

#include <cassert>
#include <vector>

namespace libint2 {

/** QuantumSet is the base class for all (sets of) quantum numbers.
    QuantumSet's must be constructable using
    std::shared_ptr<QuantumSet> or
   std::shared_ptr<ConstructablePolymorphically>.
*/
class QuantumSet : public ConstructablePolymorphically,
                   public Hashable<LIBINT2_UINT_LEAST64, ComputeKey> {
 public:
  typedef DummyIterator iter_type;
  /// Quantum numbers lie in range [0,max_quantum_number)
  static const LIBINT2_UINT_LEAST64 max_quantum_number = 100;

  virtual ~QuantumSet() {}
  virtual std::string label() const = 0;

  /// Number of quantum numbers in the set
  virtual unsigned int num_quanta() const = 0;
  /// Increment i-th quantum number
  virtual void inc(unsigned int i) = 0;
  /// Decrement i-th quantum number
  virtual void dec(unsigned int i) = 0;
};

/** QuantumNumbers<T,N> is a set of N quantum numbers of type T implemented in
 * terms of std::vector.
 */
template <typename T, unsigned int N>
class QuantumNumbers : public QuantumSet {
  std::vector<T> qn_;

 public:
  typedef QuantumSet parent_type;
  /// QuantumSet is a set of one QuantumSet
  typedef QuantumNumbers iter_type;

  QuantumNumbers(const std::vector<T>& qn);
  QuantumNumbers(const std::shared_ptr<QuantumNumbers>&);
  QuantumNumbers(const std::shared_ptr<QuantumSet>&);
  QuantumNumbers(const std::shared_ptr<ConstructablePolymorphically>&);
  QuantumNumbers(const ConstructablePolymorphically&);
  ~QuantumNumbers();

  bool operator==(const QuantumNumbers&) const;
  std::string label() const override;

  /// Increment quantum number i
  void inc(unsigned int i) override { ++qn_.at(i); }
  /// Decrement quantum number i
  void dec(unsigned int i) override {
#if CHECK_SAFETY
    if (qn_.at(i) == T(0))
      throw std::runtime_error(
          "QuantumNumbers::dec -- quantum number already zero");
#endif
    --qn_.at(i);
  }

  /// Return i-th quantum number
  const T elem(unsigned int i) const { return qn_.at(i); }

  /// Return i-th quantum number
  unsigned int num_quanta() const { return qn_.size(); }

  /// Implements Hashable::key()
  LIBINT2_UINT_LEAST64 key() const override {
    LIBINT2_UINT_LEAST64 key = 0;
    LIBINT2_UINT_LEAST64 pfac = 1;
    const int maxi = ((int)num_quanta()) - 1;
    for (int i = maxi; i >= 0; i--) {
      key += pfac * qn_[i];
      pfac *= QuantumSet::max_quantum_number;
    }
    assert(key < this->max_key());
    return key;
  }

  /// key is in range [0,max_key())
  LIBINT2_UINT_LEAST64 max_key() const {
    LIBINT2_UINT_LEAST64 max_key = 1;
    const int maxi = ((int)num_quanta()) - 1;
    for (int i = maxi; i >= 0; i--) {
      max_key *= QuantumSet::max_quantum_number;
    }
    return max_key;
  }
};

template <typename T, unsigned int N>
QuantumNumbers<T, N>::QuantumNumbers(const std::vector<T>& qn) : qn_(qn) {}

template <typename T, unsigned int N>
QuantumNumbers<T, N>::QuantumNumbers(
    const std::shared_ptr<QuantumNumbers>& sptr)
    : qn_(sptr->qn_) {}

template <typename T, unsigned int N>
QuantumNumbers<T, N>::QuantumNumbers(const std::shared_ptr<QuantumSet>& sptr) {
  const std::shared_ptr<QuantumNumbers<T, N> > sptr_cast =
      std::dynamic_pointer_cast<QuantumNumbers, QuantumSet>(sptr);
#if CHECK_SAFETY
  if (sptr_cast == 0)
    throw std::runtime_error(
        "QuantumNumbers<T,N>::QuantumNumbers(const "
        "std::shared_ptr<QuantumSet>& sptr) -- type of sptr is incompatible "
        "with QuantumNumbers");
#endif

  qn_ = sptr_cast->qn_;
}

template <typename T, unsigned int N>
QuantumNumbers<T, N>::QuantumNumbers(
    const std::shared_ptr<ConstructablePolymorphically>& sptr) {
  const std::shared_ptr<QuantumNumbers<T, N> > sptr_cast =
      std::dynamic_pointer_cast<QuantumNumbers, ConstructablePolymorphically>(
          sptr);
#if CHECK_SAFETY
  if (sptr_cast == 0)
    throw std::runtime_error(
        "QuantumNumbers<T,N>::QuantumNumbers(const "
        "std::shared_ptr<ConstructablePolymorphically>& sptr) -- type of sptr "
        "is incompatible with QuantumNumbers");
#endif

  qn_ = sptr_cast->qn_;
}

template <typename T, unsigned int N>
QuantumNumbers<T, N>::QuantumNumbers(const ConstructablePolymorphically& sptr) {
  const QuantumNumbers<T, N>& sptr_cast =
      dynamic_cast<const QuantumNumbers&>(sptr);
  qn_ = sptr_cast.qn_;
}

template <typename T, unsigned int N>
QuantumNumbers<T, N>::~QuantumNumbers() {}

template <typename T, unsigned int N>
bool QuantumNumbers<T, N>::operator==(const QuantumNumbers& a) const {
  return qn_ == a.qn_;
}

template <typename T, unsigned int N>
std::string QuantumNumbers<T, N>::label() const {
  std::ostringstream oss;
  oss << "{";
  if (qn_.size() > 0) oss << qn_[0];
  for (int i = 1; i < qn_.size(); i++) oss << "," << qn_[i];
  oss << "}";
  return oss.str();
}

// typedef QuantumNumbers<unsigned int,0> NullQuantumSet;

/**
   QuantumNumbersA<T,N> is a set of N quantum numbers of type T implemented in
   terms of a C-style array. QuantumNumbersA is faster than QuantumNumbers but
   is not as safe!
*/
template <typename T, unsigned int N>
class QuantumNumbersA : public QuantumSet {
  T qn_[N];

 public:
  typedef QuantumSet parent_type;
  /// QuantumSet is a set of one QuantumSet
  typedef QuantumNumbersA iter_type;

  // Set all quanta to val
  QuantumNumbersA(const T& val);
  QuantumNumbersA(const T* qn);
  QuantumNumbersA(const std::vector<T>& qn);
  QuantumNumbersA(const std::shared_ptr<QuantumNumbersA>&);
  QuantumNumbersA(const std::shared_ptr<QuantumSet>&);
  QuantumNumbersA(const std::shared_ptr<ConstructablePolymorphically>&);
  ~QuantumNumbersA();

  bool operator==(const QuantumNumbersA&) const;
  std::string label() const override;

  /// Increment quantum number i
  void inc(unsigned int i) override { ++qn_[i]; }
  /// Decrement quantum number i
  void dec(unsigned int i) override {
#if CHECK_SAFETY
    if (qn_[i] == T(0))
      throw std::runtime_error(
          "QuantumNumbersA::dec -- quantum number already zero");
#endif
    --qn_[i];
  }

  /// Return i-th quantum number
  const T elem(unsigned int i) const { return qn_[i]; }

  /// Return i-th quantum number
  void set_elem(unsigned int i, const T& value) { qn_[i] = value; }

  /// Implementation of QuantumSet::num_quanta()
  unsigned int num_quanta() const override { return N; }

  /// Implements Hashable::key()
  LIBINT2_UINT_LEAST64 key() const override {
    LIBINT2_UINT_LEAST64 key = 0;
    LIBINT2_UINT_LEAST64 pfac = 1;
    const int maxi = ((int)num_quanta()) - 1;
    for (int i = maxi; i >= 0; i--) {
      key += pfac * qn_[i];
      pfac *= QuantumSet::max_quantum_number;
    }
    assert(key < this->max_key());
    return key;
  }

  /// key is in range [0,max_key())
  LIBINT2_UINT_LEAST64 max_key() const {
    LIBINT2_UINT_LEAST64 max_key = 1;
    const int maxi = ((int)num_quanta()) - 1;
    for (int i = maxi; i >= 0; i--) {
      max_key *= QuantumSet::max_quantum_number;
    }
    return max_key;
  }
};

template <typename T, unsigned int N>
QuantumNumbersA<T, N>::QuantumNumbersA(const T& val) {
  for (unsigned int i = 0; i < N; i++) qn_[i] = val;
}

template <typename T, unsigned int N>
QuantumNumbersA<T, N>::QuantumNumbersA(const T* qn) {
  for (int i = 0; i < N; i++) qn_[i] = qn[i];
}

template <typename T, unsigned int N>
QuantumNumbersA<T, N>::QuantumNumbersA(const std::vector<T>& qn) {
  for (int i = 0; i < N; i++) qn_[i] = qn[i];
}

template <typename T, unsigned int N>
QuantumNumbersA<T, N>::QuantumNumbersA(
    const std::shared_ptr<QuantumNumbersA>& sptr) {
  T* qn = sptr->qn_;
  for (unsigned int i = 0; i < N; i++) qn_[i] = qn[i];
}

template <typename T, unsigned int N>
QuantumNumbersA<T, N>::QuantumNumbersA(
    const std::shared_ptr<QuantumSet>& sptr) {
  const std::shared_ptr<QuantumNumbersA<T, N> > sptr_cast =
      std::dynamic_pointer_cast<QuantumNumbersA, QuantumSet>(sptr);
#if CHECK_SAFETY
  if (sptr_cast == 0)
    throw std::runtime_error(
        "QuantumNumbersA<T,N>::QuantumNumbersA(const "
        "std::shared_ptr<QuantumSet>& sptr) -- type of sptr is incompatible "
        "with QuantumNumbersA");
#endif

  T* qn = sptr_cast->qn_;
  for (int i = 0; i < N; i++) qn_[i] = qn[i];
}

template <typename T, unsigned int N>
QuantumNumbersA<T, N>::QuantumNumbersA(
    const std::shared_ptr<ConstructablePolymorphically>& sptr) {
  const std::shared_ptr<QuantumNumbersA<T, N> > sptr_cast =
      std::dynamic_pointer_cast<QuantumNumbersA, ConstructablePolymorphically>(
          sptr);
#if CHECK_SAFETY
  if (sptr_cast == 0)
    throw std::runtime_error(
        "QuantumNumbersA<T,N>::QuantumNumbersA(const "
        "std::shared_ptr<ConstructablePolymorphically>& sptr) -- type of sptr "
        "is incompatible with QuantumNumbersA");
#endif
  T* qn = sptr_cast->qn_;
  for (int i = 0; i < N; i++) qn_[i] = qn[i];
}

template <typename T, unsigned int N>
QuantumNumbersA<T, N>::~QuantumNumbersA() {}

template <typename T, unsigned int N>
bool QuantumNumbersA<T, N>::operator==(const QuantumNumbersA& a) const {
  const T* qn0 = qn_;
  const T* qn1 = a.qn_;
  for (int i = 0; i < N; i++, ++qn0, ++qn1)
    if (*qn0 != *qn1) return false;

  return true;
}

template <typename T, unsigned int N>
std::string QuantumNumbersA<T, N>::label() const {
  std::ostringstream oss;
  oss << "{";
  if (N > 0) oss << qn_[0];
  for (unsigned int i = 1; i < N; i++) oss << "," << qn_[i];
  oss << "}";
  return oss.str();
}

/** partial specialization of QuantumNumbersA for the case N=0 */
template <typename T>
class QuantumNumbersA<T, 0> : public QuantumSet {
 public:
  typedef QuantumSet parent_type;
  /// QuantumSet is a set of one QuantumSet
  typedef QuantumNumbersA iter_type;

  QuantumNumbersA() {}
  QuantumNumbersA(const T* qn) {}
  QuantumNumbersA(const std::vector<T>& qn) {}
  QuantumNumbersA(const std::shared_ptr<QuantumNumbersA>&) {}
  QuantumNumbersA(const std::shared_ptr<QuantumSet>&) {}
  QuantumNumbersA(const std::shared_ptr<ConstructablePolymorphically>&) {}
  ~QuantumNumbersA() {}

  bool operator==(const QuantumNumbersA&) const { return true; }
  std::string label() const override { return "{}"; }

  /// Increment quantum number i
  void inc(unsigned int i) override {
    throw std::runtime_error(
        "QuantumNumbersA<T,0>::inc -- no quantum numbers to increment");
  }
  /// Decrement quantum number i
  void dec(unsigned int i) override {
    throw std::runtime_error(
        "QuantumNumbersA<T,0>::inc -- no quantum numbers to decrement");
  }
  /// Return i-th quantum number
  const T elem(unsigned int i) const {
    throw std::runtime_error(
        "QuantumNumbersA<T,0>::inc -- no quantum numbers to return");
  }
  /// Implementation of QuantumSet::num_quanta()
  unsigned int num_quanta() const override { return 0; }

  /// Implements Hashable::key()
  LIBINT2_UINT_LEAST64 key() const override { return 0; }

  /// key is in range [0,max_key())
  LIBINT2_UINT_LEAST64 max_key() const { return 1; }
};

/// Default implementation of QuantumNumbers
template <typename T, unsigned int N>
struct DefaultQuantumNumbers {
  /// This defines which QuantumNumbers implementation to use
  typedef QuantumNumbersA<T, N> Result;
};
/**
   EmptySet is the type that describes null set of auxiliary indices
*/
typedef DefaultQuantumNumbers<int, 0>::Result EmptySet;
/**
   mType is the type that describes the auxiliary index of standard 2-body
   repulsion integrals
*/
typedef DefaultQuantumNumbers<unsigned int, 1>::Result mType;

};  // namespace libint2

#endif
