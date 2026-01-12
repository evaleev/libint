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

#ifndef _libint2_src_bin_libint_algebra_h_
#define _libint2_src_bin_libint_algebra_h_

#include <class_registry.h>
#include <dgvertex.h>
#include <exception.h>
#include <global_macros.h>
#include <rr.h>
#include <smart_ptr.h>

namespace libint2 {

namespace algebra {
struct OperatorTypes {
  typedef enum { Plus, Minus, Times, Divide } OperatorType;
};
static const char OperatorSymbol[][2] = {"+", "-", "*", "/"};
};  // namespace algebra

/**
  AlgebraicOperator is an algebraic operator that
  acts on objects of type T. An algebraic operator has 2 arguments,
  to the left and to the right ( L oper R ).
*/

template <class T>
class AlgebraicOperator : public DGVertex {
 public:
  typedef algebra::OperatorTypes OperatorTypes;
  typedef algebra::OperatorTypes::OperatorType OperatorType;

  AlgebraicOperator(OperatorType OT, const std::shared_ptr<T>& left,
                    const std::shared_ptr<T>& right)
      : DGVertex(ClassInfo<AlgebraicOperator>::Instance().id()),
        OT_(OT),
        left_(left),
        right_(right),
        label_(algebra::OperatorSymbol[OT_]) {}
  virtual ~AlgebraicOperator() {}

  /// Clone A but replace operands with left and right
  AlgebraicOperator(const std::shared_ptr<AlgebraicOperator>& A,
                    const std::shared_ptr<T>& left,
                    const std::shared_ptr<T>& right)
      : DGVertex(static_cast<DGVertex&>(*A)),
        OT_(A->OT_),
        left_(left),
        right_(right),
        label_(A->label_) {
#if DEBUG
    if (num_exit_arcs() != 2)
      std::cout << "AlgebraicOperator<DGVertex> copy constructor: number of "
                   "children != 2"
                << std::endl;
    else {
      typedef DGVertex::ArcSetType::const_iterator aciter;
      aciter a = this->first_exit_arc();
      auto left_arg = (*a)->dest();
      ++a;
      auto right_arg = (*a)->dest();

      if (left_ != left_arg && left_ != right_arg)
        std::cout << "AlgebraicOperator<DGVertex> copy constructor: invalid "
                     "left operand given"
                  << std::endl;
      if (right_ != left_arg && right_ != right_arg)
        std::cout << "AlgebraicOperator<DGVertex> copy constructor: invalid "
                     "right operand given"
                  << std::endl;
    }
#endif
  }

  /// Returns the OperatorType
  OperatorType type() const { return OT_; }
  /// Returns the left argument
  const std::shared_ptr<T>& left() const { return left_; }
  /// Returns the right argument
  const std::shared_ptr<T>& right() const { return right_; }

  /// Overloads DGVertex::add_exit_arc(). The default definition is used unless
  /// T = DGVertex (see algebra.cc)
  void add_exit_arc(const std::shared_ptr<DGArc>& a) override;
  /// Implements DGVertex::size()
  unsigned int size() const override { return 1; }
  /// Implements DGVertex::equiv()
  bool equiv(const std::shared_ptr<DGVertex>& a) const override {
    if (typeid_ == a->typeid_) {
#if ALGEBRAICOPERATOR_USE_KEY_TO_COMPARE
#if USE_INT_KEY_TO_COMPARE
      if (key() == a->key())
        return *this ==
               std::static_pointer_cast<AlgebraicOperator, DGVertex>(a);
      else
        return false;
#else
      return description() == a->description();
#endif
#else
      return *this == std::static_pointer_cast<AlgebraicOperator, DGVertex>(a);
#endif
    } else
      return false;
  }

  /// laboriously compare 2 operators element by element
  bool operator==(const std::shared_ptr<AlgebraicOperator>& a) const {
#if ALGEBRAICOPERATOR_USE_SHAREDPTR
    // Find out why sometimes equivalent left_ and a->left_ have non-equivalent
    // pointers
    if (left_->equiv(a->left()) && left_ != a->left_) {
      std::cout << "Left arguments are equivalent but pointers differ!"
                << std::endl;
      std::cout << left_->description() << std::endl;
      std::cout << a->left_->description() << std::endl;
    }
    // Find out why sometimes equivalent right_ and a->right_ have
    // non-equivalent pointers
    if (right_->equiv(a->right()) && right_ != a->right_) {
      std::cout << "Left arguments are equivalent but pointers differ!"
                << std::endl;
      std::cout << right_->description() << std::endl;
      std::cout << a->right_->description() << std::endl;
    }
#endif
    if (OT_ == a->OT_) {
#if ALGEBRAICOPERATOR_USE_SHAREDPTR
      if (left_ == a->left_ && right_ == a->right_)
#else
      if (left_->equiv(a->left()) && right_->equiv(a->right()))
#endif
        return true;
      else
        return false;
    } else
      return false;
  }

  /// Implements DGVertex::label()
  const std::string& label() const override { return label_; }
  /// Implements DGVertex::id()
  const std::string& id() const override { return label(); }
  /// Implements DGVertex::description()
  std::string description() const override {
    std::ostringstream os;
    os << "( ( " << left_->description() << " ) "
       << algebra::OperatorSymbol[OT_] << " ( " << right_->description()
       << " ) )";
    const std::string descr = os.str();
    return descr;
  }

  /// Overloads DGVertex::del_exit_arcs()
  void del_exit_arcs() override {
    throw CannotAddArc(
        "AlgebraicOperator::del_exit_arcs() -- cannot safely use del_exit_arcs "
        "on operator vertices.");
  }

  /// Implements Hashable::key()
  typename DGVertex::KeyReturnType key() const override { return 0; }

  void print(std::ostream& os) const override {
    DGVertex::print(os);
    std::string prefix("AlgebraicOperator::print: ");
    os << prefix << "this = " << this << std::endl;
    os << prefix << "left_ = " << left_ << std::endl;
    os << prefix << "right_ = " << right_ << std::endl;
  }

 private:
  OperatorType OT_;
  std::shared_ptr<T> left_;
  std::shared_ptr<T> right_;

  /// Implements DGVertex::this_precomputed()
  bool this_precomputed() const override { return false; }

  std::string label_;
};

/*
template <>
  void
  AlgebraicOperator<DGVertex>::add_exit_arc(const std::shared_ptr<DGArc>& a)
  {
    DGVertex::add_exit_arc(a);
    if (left_->equiv(a->dest()))
      left_ = a->dest();
    else if (right_->equiv(a->dest()))
      right_ = a->dest();
    else
      throw std::runtime_error("AlgebraicOperator<DGVertex>::add_exit_arc --
trying to add an arc to a vertex not equivalent to either argument.");
  }
*/

/** represents linear combination of objects of type T with coefficients of type
 * C
 */
template <class C, class T>
class LinearCombination {
 public:
  typedef std::pair<C, T> term_t;
  typedef std::vector<term_t> data_t;

  LinearCombination() {}
  // shallow copy constructor -- used by operator^
  LinearCombination(data_t* data) {
    swap(*data, data_);
    delete data;
  }

  LinearCombination& operator+=(const term_t& t) {
    data_.push_back(t);
    return *this;
  }
  size_t size() const { return data_.size(); }
  const term_t& operator[](unsigned int i) const { return data_[i]; }

 private:
  data_t data_;
};

namespace algebra {
/// Wedge is a typeholder for the result of a wedge product
template <class L, class R>
struct Wedge {
  Wedge(const L& l, const R& r) : left(l), right(r) {}
  L left;
  R right;
};
template <class L, class R>
Wedge<L, R> make_wedge(const L& l, const R& r) {
  return Wedge<L, R>(l, r);
}

/// what it means to wedge terms in LinearCombination
template <class C, class Tl, class Tr>
typename LinearCombination<C, Wedge<Tl, Tr> >::term_t wedge(
    const std::pair<C, Tl>& L, const std::pair<C, Tr>& R) {
  return make_pair(L.first * R.first, L.second ^ R.second);
}
}  // namespace algebra

template <class C, class Tl, class Tr>
typename LinearCombination<C, algebra::Wedge<Tl, Tr> >::data_t* operator^(
    const LinearCombination<C, Tl>& L, const LinearCombination<C, Tr>& R) {
  typedef typename LinearCombination<C, algebra::Wedge<Tl, Tr> >::data_t data_t;
  data_t* result = new data_t;
  const size_t nL = L.size();
  const size_t nR = R.size();
  for (unsigned int l = 0; l < nL; ++l)
    for (unsigned int r = 0; r < nR; ++r) {
      result->push_back(algebra::wedge(L[l], R[r]));
    }
  return result;
}

template <class C, class Tl, class Tr>
typename LinearCombination<C, algebra::Wedge<Tl, Tr> >::data_t* operator^(
    const Tl& L, const LinearCombination<C, Tr>& R) {
  typedef typename LinearCombination<C, algebra::Wedge<Tl, Tr> >::data_t data_t;
  data_t* result = new data_t;
  const size_t nR = R.size();
  for (unsigned int r = 0; r < nR; ++r) {
    result->push_back(L ^ R[r]);
  }
  return result;
}

};  // namespace libint2

#endif
