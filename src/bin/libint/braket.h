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

#ifndef _libint2_src_bin_libint_braket_h_
#define _libint2_src_bin_libint_braket_h_

#include <algebra.h>
#include <bfset.h>
#include <global_macros.h>
#include <libint2/util/intrinsic_types.h>
#include <polyconstr.h>

namespace libint2 {

//////////////////////////////////

/** ArrayBraket is a lightweight implementation of Braket concept.
    BFS specifies the BasisFunctionSet. NP is the number of particles (>= 1).
    ArrayBraket assumes 1 BFS per particle.
*/
template <class BFS, unsigned int NP>
class ArrayBraket {
 public:
  /// This type
  typedef ArrayBraket<BFS, NP> this_type;
  /// There's no parent
  typedef struct {
  } parent_type;
  /// The iterator through ArrayBraket
  typedef ArrayBraket<typename BFS::iter_type, NP> iter_type;
  /// The total # of particles
  static constexpr auto num_particles = NP;
  /// The total # of basis functions = # of particles (since 1 bf / particle)
  static constexpr auto num_bf = NP;

  typedef BFS bfs_type;
  typedef bfs_type& bfs_ref;
  typedef const BFS& bfs_cref;

  /** This one is a very dangerous constructor -- do not to use it if at all
     possible. Provided only for compatibility for generic subiterator
     algorithms */
  ArrayBraket();
  ArrayBraket(const BFS* braket);
  ArrayBraket(const std::vector<std::vector<BFS> >& braket);
  ~ArrayBraket();

  /// Comparison function
  bool operator==(const this_type&) const;
  /// Returns cref to the i-th function for particle p
  bfs_cref member(unsigned int p, unsigned int i = 0) const;
  /// Returns ref to the i-th function for particle p
  bfs_ref member(unsigned int p, unsigned int i = 0);
  /// Returns pointer to the SubIterator for i-th BFS of particle p
  SubIterator* member_subiter(unsigned int p, unsigned int i = 0) const;
  /// Sets i-th function for particle p
  void set_member(const BFS&, unsigned int p, unsigned int i = 0);
  /// Sets i-th function for particle p (does a dynamic cast inside)
  void set_member(const ConstructablePolymorphically&, unsigned int p,
                  unsigned int i = 0);
  /// Returns the number of BFS for particle p
  unsigned int num_members(unsigned int p) const {
    assert(p < NP);
    return 1;
  }
  /// Returns the number of particles
  unsigned int num_part() const { return NP; }
  /// Implements Hashable::key()
  inline LIBINT2_UINT_LEAST64 key() const;
  /// key lies in range [0,max_key())
  LIBINT2_UINT_LEAST64 max_key() const;
#if COMPUTE_SIZE_DIRECTLY
  unsigned int size() const;
#endif

 private:
  BFS bfs_[NP];

  // this function resets all cached values
  void reset_cache();
#if COMPUTE_SIZE_DIRECTLY
  mutable unsigned int size_;
#endif
};

template <class BFS, unsigned int NP>
ArrayBraket<BFS, NP>::ArrayBraket() {
#if COMPUTE_SIZE_DIRECTLY
  size_ = 0;
#endif
}

template <class BFS, unsigned int NP>
ArrayBraket<BFS, NP>::ArrayBraket(const BFS* braket) {
  for (int i = 0; i < NP; i++) bfs_[i] = braket[i];
#if COMPUTE_SIZE_DIRECTLY
  size_ = 0;
#endif
}

template <class BFS, unsigned int NP>
ArrayBraket<BFS, NP>::ArrayBraket(
    const std::vector<std::vector<BFS> >& braket) {
  assert(braket.size() == NP);
  for (unsigned int i = 0; i < NP; i++) {
    assert(braket[i].size() == 1);
    bfs_[i] = braket[i][0];
  }
#if COMPUTE_SIZE_DIRECTLY
  size_ = 0;
#endif
}

template <class BFS, unsigned int NP>
ArrayBraket<BFS, NP>::~ArrayBraket() {}

template <class BFS, unsigned int NP>
bool ArrayBraket<BFS, NP>::operator==(const ArrayBraket<BFS, NP>& a) const {
  for (int i = 0; i < NP; i++)
    if (bfs_[i] != a[i]) return false;
  return true;
}

template <class BFS, unsigned int NP>
typename ArrayBraket<BFS, NP>::bfs_cref ArrayBraket<BFS, NP>::member(
    unsigned int p, unsigned int i) const {
  assert(i == 0 && p < NP);
  return bfs_[p];
}

template <class BFS, unsigned int NP>
typename ArrayBraket<BFS, NP>::bfs_ref ArrayBraket<BFS, NP>::member(
    unsigned int p, unsigned int i) {
  assert(i == 0 && p < NP);
  reset_cache();
  return bfs_[p];
}

template <class BFS, unsigned int NP>
SubIterator* ArrayBraket<BFS, NP>::member_subiter(unsigned int p,
                                                  unsigned int i) const {
  assert(i == 0 && p < NP);
  return static_cast<SubIterator*>(new SubIteratorBase<BFS>(member(p, i)));
}

template <class BFS, unsigned int NP>
void ArrayBraket<BFS, NP>::set_member(const BFS& f, unsigned int p,
                                      unsigned int i) {
  assert(i == 0 && p < NP);
  reset_cache();
  bfs_[p] = f;
}

template <class BFS, unsigned int NP>
void ArrayBraket<BFS, NP>::set_member(const ConstructablePolymorphically& f,
                                      unsigned int p, unsigned int i) {
  assert(i == 0 && p < NP);
  reset_cache();
  bfs_[p] = BFS(f);
}

template <class BFS, unsigned int NP>
LIBINT2_UINT_LEAST64 ArrayBraket<BFS, NP>::key() const {
  LIBINT2_UINT_LEAST64 pfac = 1;
  LIBINT2_UINT_LEAST64 key = 0;
  for (int p = NP - 1; p >= 0; p--) {
    key += pfac * bfs_[p].key();
    pfac *= BFS::max_key;
  }
  assert(key < this->max_key());
  return key;
}

template <class BFS, unsigned int NP>
LIBINT2_UINT_LEAST64 ArrayBraket<BFS, NP>::max_key() const {
  LIBINT2_UINT_LEAST64 max_key = 1;
  for (int p = NP - 1; p >= 0; p--) {
    max_key *= BFS::max_key;
  }
  return max_key;
}

#if COMPUTE_SIZE_DIRECTLY
template <class BFS, unsigned int NP>
unsigned int ArrayBraket<BFS, NP>::size() const {
  if (size_ > 0) return size_;
  size_ = 1;
  if (!TrivialBFSet<BFS>::result) {
    for (int p = NP - 1; p >= 0; p--) {
      size_ *= bfs_[p].num_bf();
    }
  }
  return size_;
}
#endif

template <class BFS, unsigned int NP>
void ArrayBraket<BFS, NP>::reset_cache() {
  size_ = 0;
}

/// This is the implementation of the Braket concept used by GenIntegralSet_1_1
// really need to have typedef template!
template <typename BFS>
struct DefaultOnePBraket {
  /// This defines which Braket implementation to use
  typedef ArrayBraket<BFS, 1> Result;
};

/// This is the implementation of the Braket concept used by
/// GenIntegralSet_11_11
// really need to have typedef template!
template <typename BFS>
struct DefaultTwoPBraket {
  /// This defines which Braket implementation to use
  typedef ArrayBraket<BFS, 2> Result;
};

///////////
/** enumerates types of brakets used for describing two-electron integrals:
    CBra and CKet are bra and ket in chemists' notation,
    PBra and PKet are bra and ket in physicists' notation.
 */
enum BraketType { CBra, CKet, PBra, PKet };
/** BraketPair is a trimmed down version of ArrayBraket specialized for
   same-particle or different-particle pairs of functions. It is not meant to be
   ArrayBraket replacement.
*/
template <class BFS, BraketType BKType>
class BraketPair {
 public:
  /// This type
  typedef BraketPair<BFS, BKType> this_type;
  typedef BFS bfs_type;

  /** This one is a very dangerous constructor -- do not to use it if at all
     possible. Provided only for compatibility for generic subiterator
     algorithms */
  BraketPair(const BFS& f1, const BFS& f2) : bfs_(f1, f2) {}
  BraketPair(const BraketPair& source) : bfs_(source.bfs_) {}
  BraketPair& operator=(const BraketPair& rhs) {
    bfs_ = rhs.bfs_;
    return *this;
  }
  const BFS& operator[](unsigned int i) const {
    if (i == 0) return bfs_.first;
    if (i == 1) return bfs_.second;
    throw std::logic_error("BraketPair::operator[] -- argument out of range");
  }
  /// Comparison function
  bool operator==(const this_type& rhs) const { return rhs.bfs_ == bfs_; }

 private:
  std::pair<BFS, BFS> bfs_;
};

/// these objects help to construct BraketPairs
namespace braket {
/// Physicists bra
template <class F>
BraketPair<F, PBra> _pbra(const F& f1, const F& f2) {
  return BraketPair<F, PBra>(f1, f2);
}
/// Physicists ket
template <class F>
BraketPair<F, PKet> _pket(const F& f1, const F& f2) {
  return BraketPair<F, PKet>(f1, f2);
}
/// Chemists bra
template <class F>
BraketPair<F, CBra> _cbra(const F& f1, const F& f2) {
  return BraketPair<F, CBra>(f1, f2);
}
/// Chemists ket
template <class F>
BraketPair<F, CKet> _cket(const F& f1, const F& f2) {
  return BraketPair<F, CKet>(f1, f2);
}
};  // namespace braket

template <class BFS, BraketType BKTypeL, BraketType BKTypeR>
algebra::Wedge<BraketPair<BFS, BKTypeL>, BraketPair<BFS, BKTypeR> > operator^(
    const BraketPair<BFS, BKTypeL>& L, const BraketPair<BFS, BKTypeR>& R) {
  return algebra::make_wedge(L, R);
}

};  // namespace libint2

#endif
