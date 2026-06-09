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

#ifndef _libint2_src_bin_libint_compderivgaussv2_h_
#define _libint2_src_bin_libint_compderivgaussv2_h_

#include <algebra.h>
#include <context.h>
#include <default_params.h>
#include <dgvertex.h>
#include <dims.h>
#include <integral.h>
#include <prefactors.h>
#include <rr.h>
#include <task.h>

#include <cassert>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace libint2 {

/** Optimized compute relation for (geometric) derivative Gaussian integrals.
 *
 * Like CR_DerivGauss, this expands derivative Gaussians via:
 *   d/dr G(a) = 2*alpha * G(a+1) - a_i * G(a-1)
 *
 * Unlike CR_DerivGauss, this uses the HRR-like code-sharing optimization:
 * since differentiation of a Gaussian at position (part, where) depends only
 * on that shell's quanta (not spectator shells), we generate code once per
 * unique differentiated shell and pass spectator dimensions at the call site.
 *
 * @tparam IntType integral type
 * @tparam part particle index of the function to be differentiated
 * @tparam where position of the function to be differentiated (InBra/InKet)
 * @tparam trans_inv_part if non-negative, specifies the particle index for
 *         translational invariance
 * @tparam trans_inv_where position for translational invariance
 */
template <class IntType, int part, FunctionPosition where,
          int trans_inv_part = -1, FunctionPosition trans_inv_where = InBra>
class DerivGaussV2 : public RecurrenceRelation {
 private:
  static constexpr auto trans_inv_oper =
      not IntType::OperType::Properties::odep;
  static constexpr auto using_trans_inv =
      trans_inv_oper && (part == trans_inv_part) && (where == trans_inv_where);

 public:
  typedef RecurrenceRelation ParentType;
  typedef typename IntType::BasisFunctionType BasisFunctionType;
  typedef DerivGaussV2 ThisType;
  typedef IntType TargetType;
  typedef IntType ChildType;
  typedef RecurrenceRelation::ExprType ExprType;

  static const unsigned int max_nchildren_ =
      using_trans_inv ? (IntType::num_bf - 1) : 2u;

  static std::shared_ptr<ThisType> Instance(
      const std::shared_ptr<TargetType>& Tint, unsigned int dir = 0);
  virtual ~DerivGaussV2() {}

  /// always directional
  static bool directional() { return true; }

  unsigned int num_children() const override { return nchildren_; }
  std::shared_ptr<DGVertex> rr_target() const override {
    return std::static_pointer_cast<DGVertex, TargetType>(target_);
  }
  std::shared_ptr<DGVertex> rr_child(unsigned int i) const override {
    return children_.at(i);
  }
  bool is_simple() const override {
    return TrivialBFSet<BasisFunctionType>::result;
  }

  std::string spfunction_call(
      const std::shared_ptr<CodeContext>& context,
      const std::shared_ptr<ImplicitDimensions>& dims) const override;

 private:
  DerivGaussV2(const std::shared_ptr<TargetType>& Tint, unsigned int dir);

  unsigned int dir_;
  std::shared_ptr<TargetType> target_;
  std::vector<std::shared_ptr<DGVertex>> children_;
  unsigned int nchildren_;

  std::string generate_label() const override;
  std::shared_ptr<ImplicitDimensions> adapt_dims_(
      const std::shared_ptr<ImplicitDimensions>& dims) const override;
  bool register_with_rrstack() const;
  bool expl_high_dim() const;
  bool expl_low_dim() const;

  /// add child, deduplicating
  const std::shared_ptr<DGVertex>& add_child(
      const std::shared_ptr<DGVertex>& child) {
    for (auto& c : children_) {
      if (c == child) return c;
    }
    children_.push_back(child);
    ++nchildren_;
    return children_.back();
  }
};

//
// Implementation
//

template <class IntType, int part, FunctionPosition where, int trans_inv_part,
          FunctionPosition trans_inv_where>
std::shared_ptr<
    DerivGaussV2<IntType, part, where, trans_inv_part, trans_inv_where>>
DerivGaussV2<IntType, part, where, trans_inv_part, trans_inv_where>::Instance(
    const std::shared_ptr<TargetType>& Tint, unsigned int dir) {
  std::shared_ptr<ThisType> this_ptr(new ThisType(Tint, dir));
  if (this_ptr->num_children() != 0) {
    this_ptr->register_with_rrstack();
    return this_ptr;
  }
  return std::shared_ptr<ThisType>();
}

template <class IntType, int part, FunctionPosition where, int trans_inv_part,
          FunctionPosition trans_inv_where>
DerivGaussV2<IntType, part, where, trans_inv_part, trans_inv_where>::
    DerivGaussV2(const std::shared_ptr<TargetType>& Tint, unsigned int dir)
    : dir_(dir), target_(Tint), nchildren_(0) {
  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using namespace libint2::braket;
  typedef BasisFunctionType F;
  const F& _1 = unit<F>(is_simple() ? dir : 0);

  const typename IntType::AuxQuantaType& aux = Tint->aux();
  const typename IntType::OperType& oper = Tint->oper();

  children_.reserve(max_nchildren_);

  // the Gaussian must be differentiated in direction dir
  {
    if (where == InBra && Tint->bra(part, 0).deriv().d(dir) == 0) return;
    if (where == InKet && Tint->ket(part, 0).deriv().d(dir) == 0) return;
  }

  // if not using translational invariance, can only expand primitives
  if (not using_trans_inv) {
    if (where == InBra && Tint->bra(part, 0).contracted()) return;
    if (where == InKet && Tint->ket(part, 0).contracted()) return;
  }

  typedef typename IntType::BraType IBraType;
  typedef typename IntType::KetType IKetType;
  IBraType* bra = new IBraType(Tint->bra());
  IKetType* ket = new IKetType(Tint->ket());

  if (not using_trans_inv) {  // differentiate

    if (where == InBra) {
      F a(bra->member(part, 0));

      // add a+1
      F ap1(bra->member(part, 0) + _1);
      ap1.deriv().dec(dir);
      bra->set_member(ap1, part, 0);
      auto int_ap1 = add_child(IntType::Instance(*bra, *ket, aux, oper));
      bra->set_member(a, part, 0);
      if (is_simple()) {
        std::ostringstream oss;
        oss << "two_alpha" << part << "_bra";
        expr_ = Scalar(oss.str()) * int_ap1;
        nflops_ += 1;
      }

      // See if a-1 exists
      F am1(bra->member(part, 0) - _1);
      if (exists(am1)) {
        am1.deriv().dec(dir);
        bra->set_member(am1, part, 0);
        auto int_am1 = add_child(IntType::Instance(*bra, *ket, aux, oper));
        bra->set_member(a, part, 0);
        if (is_simple()) {
          expr_ -= Scalar(a[dir]) * int_am1;
          nflops_ += 2;
        }
      }
      delete bra;
      delete ket;
      return;
    }

    if (where == InKet) {
      F a(ket->member(part, 0));

      // add a+1
      F ap1(ket->member(part, 0) + _1);
      ap1.deriv().dec(dir);
      ket->set_member(ap1, part, 0);
      auto int_ap1 = add_child(IntType::Instance(*bra, *ket, aux, oper));
      ket->set_member(a, part, 0);
      if (is_simple()) {
        std::ostringstream oss;
        oss << "two_alpha" << part << "_ket";
        expr_ = Scalar(oss.str()) * int_ap1;
        nflops_ += 1;
      }

      // See if a-1 exists
      F am1(ket->member(part, 0) - _1);
      if (exists(am1)) {
        am1.deriv().dec(dir);
        ket->set_member(am1, part, 0);
        auto int_am1 = add_child(IntType::Instance(*bra, *ket, aux, oper));
        ket->set_member(a, part, 0);
        if (is_simple()) {
          expr_ -= Scalar(a[dir]) * int_am1;
          nflops_ += 2;
        }
      }
      delete bra;
      delete ket;
      return;
    }

  } else {  // use translational invariance

    // remove one deriv quantum from the target function
    if (where == InBra) bra->member(part, 0).deriv().dec(dir);
    if (where == InKet) ket->member(part, 0).deriv().dec(dir);

    int term_count = 0;
    for (int p = 0; p != IntType::num_particles; ++p) {
      typedef BasisFunctionType F;
      if (p != trans_inv_part || trans_inv_where != InBra) {
        F a(bra->member(p, 0));
        if (not a.is_unit()) {
          F da(a);
          da.deriv().inc(dir);
          bra->set_member(da, p, 0);
          auto int_da = add_child(IntType::Instance(*bra, *ket, aux, oper));
          bra->set_member(a, p, 0);
          if (is_simple()) {
            if (term_count == 0)
              expr_ = Scalar(-1) * int_da;
            else
              expr_ -= int_da;
            ++term_count;
            nflops_ += 1;
          }
        }
      }
      if (p != trans_inv_part || trans_inv_where != InKet) {
        F a(ket->member(p, 0));
        if (not a.is_unit()) {
          F da(a);
          da.deriv().inc(dir);
          ket->set_member(da, p, 0);
          auto int_da = add_child(IntType::Instance(*bra, *ket, aux, oper));
          ket->set_member(a, p, 0);
          if (is_simple()) {
            if (term_count == 0)
              expr_ = Scalar(-1) * int_da;
            else
              expr_ -= int_da;
            ++term_count;
            nflops_ += 1;
          }
        }
      }
    }
  }

  delete bra;
  delete ket;
}

template <class IntType, int part, FunctionPosition where, int trans_inv_part,
          FunctionPosition trans_inv_where>
bool DerivGaussV2<IntType, part, where, trans_inv_part,
                  trans_inv_where>::register_with_rrstack() const {
  using std::swap;

  // only register RRs for shell sets (not individual integrals)
  if (TrivialBFSet<BasisFunctionType>::result) return false;

  // translational invariance path not optimized yet — register as-is
  if (using_trans_inv) {
    std::shared_ptr<RRStack> rrstack = RRStack::Instance();
    std::shared_ptr<ThisType> this_ptr =
        std::const_pointer_cast<ThisType, const ThisType>(
            std::static_pointer_cast<const ThisType, const ParentType>(
                std::enable_shared_from_this<ParentType>::shared_from_this()));
    rrstack->find(this_ptr);
    return true;
  }

  typedef typename IntType::BraType IBraType;
  typedef typename IntType::KetType IKetType;
  const IBraType& bra = target_->bra();
  const IKetType& ket = target_->ket();

  // check if all spectator shells already have zero quanta
  bool nonzero_quanta = false;
  unsigned const int npart = IntType::OperatorType::Properties::np;
  for (unsigned int p = 0; p < npart; p++) {
    int nfbra = bra.num_members(p);
    for (int f = 0; f < nfbra; f++) {
      // skip the differentiated position
      if (static_cast<int>(p) == part && where == InBra) continue;
      if (!bra.member(p, f).zero() || !bra.member(p, f).deriv().zero())
        nonzero_quanta = true;
    }
    int nfket = ket.num_members(p);
    for (int f = 0; f < nfket; f++) {
      if (static_cast<int>(p) == part && where == InKet) continue;
      if (!ket.member(p, f).zero() || !ket.member(p, f).deriv().zero())
        nonzero_quanta = true;
    }
  }

  // if all spectators are zero, register this instance directly
  if (!nonzero_quanta) {
    std::shared_ptr<RRStack> rrstack = RRStack::Instance();
    std::shared_ptr<ThisType> this_ptr =
        std::const_pointer_cast<ThisType, const ThisType>(
            std::static_pointer_cast<const ThisType, const ParentType>(
                std::enable_shared_from_this<ParentType>::shared_from_this()));
    rrstack->find(this_ptr);
    return true;
  }

  // Otherwise, zero out all spectator shells and register a dummy
  IBraType bra_zero(bra);
  IKetType ket_zero(ket);
  for (unsigned int p = 0; p < npart; p++) {
    int nfbra = bra_zero.num_members(p);
    for (int f = 0; f < nfbra; f++) {
      if (static_cast<int>(p) == part && where == InBra) continue;
      typedef typename IBraType::bfs_type bfs_type;
      typedef typename IBraType::bfs_ref bfs_ref;
      bfs_ref bfs = bra_zero.member(p, f);
      if (!bfs.zero() || !bfs.deriv().zero()) {
        bfs_type null_bfs;
        swap(bfs, null_bfs);
      }
    }
    int nfket = ket_zero.num_members(p);
    for (int f = 0; f < nfket; f++) {
      if (static_cast<int>(p) == part && where == InKet) continue;
      typedef typename IKetType::bfs_type bfs_type;
      typedef typename IKetType::bfs_ref bfs_ref;
      bfs_ref bfs = ket_zero.member(p, f);
      if (!bfs.zero() || !bfs.deriv().zero()) {
        bfs_type null_bfs;
        swap(bfs, null_bfs);
      }
    }
  }

  // create a generic integral with a dummy operator
  typedef GenOper<GenMultSymmOper_Descr<IntType::OperatorType::Properties::np>>
      DummyOper;
  typedef EmptySet DummyQuanta;
  typedef GenIntegralSet<DummyOper, IncableBFSet, IBraType, IKetType,
                         DummyQuanta>
      DummyIntegral;
  DummyOper dummy_oper;
  DummyQuanta dummy_quanta(std::vector<int>(0, 0));
  std::shared_ptr<DummyIntegral> dummy_integral =
      DummyIntegral::Instance(bra_zero, ket_zero, dummy_quanta, dummy_oper);

  // construct a DerivGaussV2 over the dummy integral and register it
  typedef DerivGaussV2<DummyIntegral, part, where> DummyDerivGaussV2;
  std::shared_ptr<DummyDerivGaussV2> dummy_rr =
      DummyDerivGaussV2::Instance(dummy_integral, dir_);
  std::shared_ptr<RRStack> rrstack = RRStack::Instance();
  rrstack->find(dummy_rr);
  return true;
}

template <class IntType, int part, FunctionPosition where, int trans_inv_part,
          FunctionPosition trans_inv_where>
std::string DerivGaussV2<IntType, part, where, trans_inv_part,
                         trans_inv_where>::generate_label() const {
  std::ostringstream os;

  // For translational invariance, children depend on ALL shells, so
  // the label must include full integral info (no code sharing).
  // For direct differentiation, only the differentiated shell matters.
  if constexpr (using_trans_inv) {
    typedef typename TargetType::AuxIndexType mType;
    static std::shared_ptr<mType> aux0(new mType(0u));
    os << "CR_DerivGauss"
       << "P" << part << to_string(where)
       << genintegralset_label(target_->bra(), target_->ket(), aux0,
                               target_->oper());
    return os.str();
  }

  os << "DerivGaussV2 P" << part << " " << to_string(where) << " ";

  // Only encode the differentiated shell — not spectators
  if (where == InBra) {
    BasisFunctionType sh(target_->bra(part, 0));
    sh.uncontract();
    os << sh.label();
  } else {
    BasisFunctionType sh(target_->ket(part, 0));
    sh.uncontract();
    os << sh.label();
  }

  return os.str();
}

template <class IntType, int part, FunctionPosition where, int trans_inv_part,
          FunctionPosition trans_inv_where>
std::string
DerivGaussV2<IntType, part, where, trans_inv_part, trans_inv_where>::
    spfunction_call(const std::shared_ptr<CodeContext>& context,
                    const std::shared_ptr<ImplicitDimensions>& dims) const {
  std::ostringstream os;
  os << context->label_to_function_name(label()) << "(inteval, "
     << context->value_to_pointer(rr_target()->symbol());

  const unsigned int nc = num_children();
  for (unsigned int c = 0; c < nc; c++) {
    os << ", " << context->value_to_pointer(rr_child(c)->symbol());
  }

  // compute hsr and lsr — dimensions of spectator shells
  // canonical order: for each particle p, bra then ket
  // hsr = product of dims before (part, where)
  // lsr = product of dims after (part, where)
  unsigned int hsr = 1;
  unsigned int lsr = 1;
  const unsigned int np = IntType::OperType::Properties::np;
  for (int p = 0; p < static_cast<int>(np); p++) {
    unsigned int nbra = target_->bra().num_members(p);
    assert(nbra == 1);
    for (unsigned int i = 0; i < nbra; i++) {
      SubIterator* iter = target_->bra().member_subiter(p, i);
      if (p < part || (p == part && where == InKet)) hsr *= iter->num_iter();
      // skip p == part && where == InBra (the differentiated shell)
      if (p > part) lsr *= iter->num_iter();
      delete iter;
    }
    unsigned int nket = target_->ket().num_members(p);
    assert(nket == 1);
    for (unsigned int i = 0; i < nket; i++) {
      SubIterator* iter = target_->ket().member_subiter(p, i);
      if (p < part) hsr *= iter->num_iter();
      // skip p == part && where == InKet (the differentiated shell)
      if (p > part || (p == part && where == InBra)) lsr *= iter->num_iter();
      delete iter;
    }
  }

  // Use TaskParameters to keep track of maximum ranks
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current().params()->max_hrr_hsrank(hsr);

  if (expl_high_dim()) os << "," << hsr;
  if (expl_low_dim()) os << "," << lsr;
  os << ")" << context->end_of_stat() << std::endl;
  return os.str();
}

template <class IntType, int part, FunctionPosition where, int trans_inv_part,
          FunctionPosition trans_inv_where>
bool DerivGaussV2<IntType, part, where, trans_inv_part,
                  trans_inv_where>::expl_high_dim() const {
  // translational invariance: no code sharing, no explicit dims
  if (using_trans_inv) return false;
  // need explicit high dim unless this is the first position
  if (part == 0 && where == InBra) return false;
  return true;
}

template <class IntType, int part, FunctionPosition where, int trans_inv_part,
          FunctionPosition trans_inv_where>
bool DerivGaussV2<IntType, part, where, trans_inv_part,
                  trans_inv_where>::expl_low_dim() const {
  // translational invariance: no code sharing, no explicit dims
  if (using_trans_inv) return false;
  unsigned int np = IntType::OperType::Properties::np;
  // need explicit low dim unless this is the last position
  if (static_cast<int>(np) - 1 == part && where == InKet) return false;
  // corner case: 1-particle operator
  if (np == 1) return true;
  return true;
}

template <class IntType, int part, FunctionPosition where, int trans_inv_part,
          FunctionPosition trans_inv_where>
std::shared_ptr<ImplicitDimensions>
DerivGaussV2<IntType, part, where, trans_inv_part, trans_inv_where>::
    adapt_dims_(const std::shared_ptr<ImplicitDimensions>& dims) const {
  bool high_rank = expl_high_dim();
  bool low_rank = expl_low_dim();

  std::shared_ptr<Entity> high_dim, low_dim;
  if (high_rank) {
    high_dim =
        std::shared_ptr<Entity>(new RTimeEntity<EntityTypes::Int>("highdim"));
  } else {
    high_dim = dims->high();
  }
  if (low_rank) {
    low_dim =
        std::shared_ptr<Entity>(new RTimeEntity<EntityTypes::Int>("lowdim"));
  } else {
    low_dim = dims->low();
  }

  std::shared_ptr<ImplicitDimensions> localdims(
      new ImplicitDimensions(high_dim, low_dim, dims->vecdim()));
  return localdims;
}

};  // namespace libint2

#endif
