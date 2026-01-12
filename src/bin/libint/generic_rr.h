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

#ifndef _libint2_src_bin_libint_genericrr_h_
#define _libint2_src_bin_libint_genericrr_h_

#include <algebra.h>
#include <context.h>
#include <default_params.h>
#include <dgvertex.h>
#include <flop.h>
#include <integral.h>
#include <prefactors.h>
#include <rr.h>
#include <util.h>

#include <boost/type_traits/is_same.hpp>
#include <cassert>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace libint2 {

/** RRImpl must inherit GenericRecurrenceRelation<RRImpl>
 */
template <typename RRImpl, typename F, typename Target>
class GenericRecurrenceRelation : public RecurrenceRelation {
 public:
  typedef F BasisFunctionType;
  typedef Target TargetType;
  typedef RecurrenceRelation ParentType;
  typedef ParentType::ExprType ExprType;

  /// Return an instance if applicable, or a null pointer otherwise
  static std::shared_ptr<RRImpl> Instance(
      const std::shared_ptr<TargetType>& Tint, unsigned int dir) {
    // screen out calls with nondefault extra parameters
    if (!RRImpl::directional() && dir != 0) return std::shared_ptr<RRImpl>();
    // attempt to construct
    std::shared_ptr<RRImpl> this_ptr(new RRImpl(Tint, dir));
    // if succeeded (nchildren > 0) do post-construction
    if (this_ptr->num_children() != 0) {
      this_ptr->template register_with_rrstack<RRImpl>();
      return this_ptr;
    }
    // else return null pointer
    return std::shared_ptr<RRImpl>();
  }

  /// Implementation of RecurrenceRelation::num_children()
  unsigned int num_children() const override { return children_.size(); }
  /// Implementation of RecurrenceRelation::rr_target()
  std::shared_ptr<DGVertex> rr_target() const override {
    return std::static_pointer_cast<DGVertex, TargetType>(target_);
  }
  /// Implementation of RecurrenceRelation::rr_child()
  std::shared_ptr<DGVertex> rr_child(unsigned int i) const override {
    return children_.at(i);
  }
  /// Implementation of RecurrenceRelation::is_simple()
  bool is_simple() const override {
    return TrivialBFSet<BasisFunctionType>::result;
  }

  /// Implementation of RecurrenceRelation::generate_label()
  std::string generate_label() const override {
    std::ostringstream os;
    os << RRImpl::descr() << " " << target_->label();
    return os.str();
  }

 protected:
  GenericRecurrenceRelation(const std::shared_ptr<TargetType>& Tint,
                            unsigned int dir)
      : target_(Tint), dir_(dir) {
    children_.reserve(RRImpl::max_nchildren);
  }
  /** is this recurrence relation parameterized by a direction (x, y, or z).
      true if BasisFunctionSet is CGF,
      false otherwise. */
  static bool default_directional() {
    if (boost::is_same<BasisFunctionType, CGF>::value) return true;
    return false;
  }

  unsigned int dir() const { return dir_; }

  /// add child
  const std::shared_ptr<DGVertex>& add_child(
      const std::shared_ptr<DGVertex>& child) {
    typedef std::vector<std::shared_ptr<DGVertex> > cvector;
    typedef typename cvector::const_iterator citer;
    const citer pos = std::find(children_.begin(), children_.end(), child);
    if (pos == children_.end()) {
      children_.push_back(child);
      return *(children_.rbegin());
    } else
      return *pos;
  }

  /// use this helper to make children
  template <class RR, class C>
  friend class ChildFactory;
  /// make_child should really looks something like this, but gcc 4.3.0 craps
  /// out
  /// TODO test is this works
#if 0
      template <class RealChildType>
      const std::shared_ptr<DGVertex>& make_child(const typename RealChildType::BasisFunctionType& A,
                                          const typename RealChildType::BasisFunctionType& B,
                                          const typename RealChildType::BasisFunctionType& C,
                                          const typename RealChildType::BasisFunctionType& D,
                                          const typename RealChildType::AuxIndexType& aux = typename RealChildType::AuxIndexType(),
                                          const typename RealChildType::OperType& oper = typename RealChildType::OperType()) {
        const std::shared_ptr<DGVertex>& i = std::static_pointer_cast<DGVertex,RealChildType>(ChildType::Instance(A,B,C,D,aux,oper));
        return add_child(i);
      }
#endif

  std::shared_ptr<TargetType> target_;

 private:
  unsigned int dir_;
  std::vector<std::shared_ptr<DGVertex> > children_;
};

/// Helps GenericRecurrenceRelation to work around the compiler problem with
/// make_child
template <class GenRR, class ChildType>
class ChildFactory {
 public:
  typedef typename ChildType::BasisFunctionType F;
  typedef typename ChildType::AuxIndexType AuxIndexType;
  typedef typename ChildType::OperType OperType;

  ChildFactory(GenRR* rr) : rr_(rr) {}

  /// make_child
  const std::shared_ptr<DGVertex>& make_child(
      const F& A, const F& B, const F& C, const F& D,
      const AuxIndexType& aux = AuxIndexType(),
      const OperType& oper = OperType()) {
    auto i = std::static_pointer_cast<DGVertex, ChildType>(
        ChildType::Instance(A, B, C, D, aux, oper));
    return rr_->add_child(i);
  }
  /// make_child
  const std::shared_ptr<DGVertex>& make_child(
      const F& A, const F& B, const AuxIndexType& aux = AuxIndexType(),
      const OperType& oper = OperType()) {
    auto i = std::static_pointer_cast<DGVertex, ChildType>(
        ChildType::Instance(A, B, aux, oper));
    return rr_->add_child(i);
  }
  /// make a child from a wedge of physicists' brackets
  const std::shared_ptr<DGVertex>& make_child(
      const algebra::Wedge<BraketPair<F, PBra>, BraketPair<F, PKet> >&
          braket_wedge,
      const AuxIndexType& aux = AuxIndexType(),
      const OperType& oper = OperType()) {
    auto i = std::static_pointer_cast<DGVertex, ChildType>(
        ChildType::Instance(braket_wedge, aux, oper));
    return rr_->add_child(i);
  }
  /// make a child from a wedge of chemists' brackets
  const std::shared_ptr<DGVertex>& make_child(
      const algebra::Wedge<BraketPair<F, CBra>, BraketPair<F, CKet> >&
          braket_wedge,
      const AuxIndexType& aux = AuxIndexType(),
      const OperType& oper = OperType()) {
    auto i = std::static_pointer_cast<DGVertex, ChildType>(
        ChildType::Instance(braket_wedge, aux, oper));
    return rr_->add_child(i);
  }
  /// take a wedge product of various (linear combinations of) brakets
  void wedge(const LinearCombination<std::shared_ptr<DGVertex>,
                                     BraketPair<F, PBra> >& bra_lc,
             const LinearCombination<std::shared_ptr<DGVertex>,
                                     BraketPair<F, PKet> >& ket_lc,
             const AuxIndexType& aux = AuxIndexType(),
             const OperType& oper = OperType()) {
    using namespace libint2::algebra;
    typedef LinearCombination<std::shared_ptr<DGVertex>,
                              Wedge<BraketPair<F, PBra>, BraketPair<F, PKet> > >
        ProductLC;
    const ProductLC& product_lc = bra_lc ^ ket_lc;
    const size_t nprod = product_lc.size();
    for (unsigned int t = 0; t < nprod; ++t) {
      const typename ProductLC::term_t& term = product_lc[t];
      auto child = make_child(term.second, aux, oper);
      if (rr_->is_simple()) {
        if (rr_->expr_)
          rr_->expr_ += term.first * child;
        else
          rr_->expr_ = term.first * child;
      }
    }
  }
  void wedge(const BraketPair<F, PBra>& bra,
             const LinearCombination<std::shared_ptr<DGVertex>,
                                     BraketPair<F, PKet> >& ket_lc,
             const AuxIndexType& aux = AuxIndexType(),
             const OperType& oper = OperType()) {
    using namespace libint2::prefactor;
    LinearCombination<std::shared_ptr<DGVertex>, BraketPair<F, PBra> > bra_lc;
    bra_lc += make_pair(Scalar(1.0), bra);
    wedge(bra_lc, ket_lc, aux, oper);
  }
  void wedge(const LinearCombination<std::shared_ptr<DGVertex>,
                                     BraketPair<F, PBra> >& bra_lc,
             const BraketPair<F, PKet>& ket,
             const AuxIndexType& aux = AuxIndexType(),
             const OperType& oper = OperType()) {
    using namespace libint2::prefactor;
    LinearCombination<std::shared_ptr<DGVertex>, BraketPair<F, PKet> > ket_lc;
    ket_lc += make_pair(Scalar(1.0), ket);
    wedge(bra_lc, ket_lc, aux, oper);
  }

 private:
  GenRR* rr_;
};

};  // namespace libint2

#endif
