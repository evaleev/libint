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

#ifndef _libint2_src_bin_libint_hrr_h_
#define _libint2_src_bin_libint_hrr_h_

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
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace libint2 {

/** A generic Horizontal Recurrence Relation:

 |a b) = |a+1_i b-1_i) + AB_i |a b-1_i) + extra terms for derivative Gaussians

 Int is the integral class. part specifies for which particle
 the angular momentum is shifted. Function a is assumed to gain quanta,
 function b loses quanta. loc_a and loc_b specify where
 functions a and b are located (bra or ket). pos_a and pos_b
 specify which function to be used (usually pos_a and pos_b are set
 to 0 to refer to the first function for this particle in this location).

 */
template <class IntType, class BFSet, int part, FunctionPosition loc_a,
          unsigned int pos_a, FunctionPosition loc_b, unsigned int pos_b>
class HRR : public RecurrenceRelation {
 public:
  typedef RecurrenceRelation ParentType;
  typedef BFSet BasisFunctionType;
  typedef HRR<IntType, BFSet, part, loc_a, pos_a, loc_b, pos_b> ThisType;
  typedef IntType TargetType;
  typedef IntType ChildType;
  /// A short alias
  typedef RecurrenceRelation::ExprType ExprType;

  /** Use Instance() to obtain an instance of RR. This function is provided to
   avoid issues with getting a std::shared_ptr from constructor (as needed for
   registry to work).

   dir specifies which quantum number of a and b is shifted.
   For example, dir can be 0 (x), 1(y), or 2(z) if F is
   a Cartesian Gaussian.
   */
  static std::shared_ptr<ThisType> Instance(const std::shared_ptr<TargetType>&,
                                            unsigned int dir = 0);
  virtual ~HRR();

  /// overrides RecurrenceRelation::braket_direction()
  BraketDirection braket_direction() const override {
    if (loc_b == InBra && loc_a == InKet)
      return BraketDirection::BraToKet;
    else if (loc_b == InKet && loc_a == InBra)
      return BraketDirection::KetToBra;
    else
      return BraketDirection::None;
  }

  /** is this recurrence relation parameterized by a direction.
      the default is false if BasisFunctionSet is CGShell,
      true otherwise. */
  static bool directional() {
    if (boost::is_same<BasisFunctionType, CGShell>::value) return false;
    return true;
  }

  /// Implementation of RecurrenceRelation::num_children()
  unsigned int num_children() const override { return nchildren_; }
  /// returns pointer to the target
  std::shared_ptr<TargetType> target() const { return target_; };
  /// child(i) returns pointer to i-th child
  std::shared_ptr<ChildType> child(unsigned int i) const;
  /// Implementation of RecurrenceRelation::target()
  std::shared_ptr<DGVertex> rr_target() const override {
    return std::static_pointer_cast<DGVertex, TargetType>(target());
  }
  /// Implementation of RecurrenceRelation::child()
  std::shared_ptr<DGVertex> rr_child(unsigned int i) const override {
    return std::static_pointer_cast<DGVertex, ChildType>(child(i));
  }
  /// Implementation of RecurrenceRelation::is_simple()
  bool is_simple() const override { return TrivialBFSet<BFSet>::result; }
  /// Implementation of RecurrenceRelation::spfunction_call()
  std::string spfunction_call(
      const std::shared_ptr<CodeContext>& context,
      const std::shared_ptr<ImplicitDimensions>& dims) const override;

 private:
  /**
   dir specifies which quantum number of a and b is shifted.
   For example, dir can be 0 (x), 1(y), or 2(z) if F is
   a Cartesian Gaussian.
   */
  HRR(const std::shared_ptr<TargetType>&, unsigned int dir);

  unsigned int dir_;
  std::shared_ptr<TargetType> target_;
  static const unsigned int max_nchildren_ = 8;
  std::shared_ptr<ChildType> children_[max_nchildren_];
  unsigned int nchildren_;

  void oper_checks() const;

  /// Implementation of RecurrenceRelation::label()
  std::string generate_label() const override;
  /// Reimplementation of RecurrenceRelation::adapt_dims_()
  std::shared_ptr<ImplicitDimensions> adapt_dims_(
      const std::shared_ptr<ImplicitDimensions>& dims) const override;
  /// Use instead of RecurrenceRelation::register_with_rrstack()
  bool register_with_rrstack() const;
  /** return true if the high dimension must be shown explicitly. For example,
   cd-HRR applied (ss|pp) has high dimension of rank 1 but since the code for
   such RR is not specific to ab=(ss|, the rank of high dimension must be shown
   explicitly.
   */
  bool expl_high_dim() const;
  bool expl_low_dim() const;
};

template <class IntType, class F, int part, FunctionPosition loc_a,
          unsigned int pos_a, FunctionPosition loc_b, unsigned int pos_b>
std::shared_ptr<HRR<IntType, F, part, loc_a, pos_a, loc_b, pos_b> >
HRR<IntType, F, part, loc_a, pos_a, loc_b, pos_b>::Instance(
    const std::shared_ptr<TargetType>& Tint, unsigned int dir) {
  std::shared_ptr<ThisType> this_ptr(new ThisType(Tint, dir));
  // Do post-construction duties
  if (this_ptr->num_children() != 0) {
    this_ptr->register_with_rrstack();
    return this_ptr;
  }
  return std::shared_ptr<ThisType>();
}

template <class IntType, class F, int part, FunctionPosition loc_a,
          unsigned int pos_a, FunctionPosition loc_b, unsigned int pos_b>
HRR<IntType, F, part, loc_a, pos_a, loc_b, pos_b>::HRR(
    const std::shared_ptr<TargetType>& Tint, unsigned int dir)
    : dir_(dir), target_(Tint), nchildren_(0) {
  using namespace libint2::algebra;
  using namespace libint2::prefactor;

  // assume that always transfering from Bra to Ket and vice versa
  assert(loc_a != loc_b);

  target_ = Tint;
  const typename IntType::AuxQuantaType& aux = Tint->aux();
  const typename IntType::OperType& oper = Tint->oper();

  // can move across operator only if it's multiplicative
  if (loc_a != loc_b && oper.hermitian(part) != +1) {
    return;
  }

  typedef typename IntType::BraType IBraType;
  typedef typename IntType::KetType IKetType;
  IBraType* bra = new IBraType(Tint->bra());
  IKetType* ket = new IKetType(Tint->ket());

  //
  // InBra and InKet cases have to be treated explicitly since BraType and
  // KetType don't have to match
  //
  if (loc_a == InKet && loc_b == InBra) {
    F a(ket->member(part, pos_a));
    F b(bra->member(part, pos_b));

    // See if b-1 exists
    F bm1(b);
    bm1.dec(dir_);
    if (!exists(bm1)) {
      delete bra;
      delete ket;
      return;
    }
    bra->set_member(bm1, part, pos_b);  // set b permanently to b-1_i

    {
      F ap1(a);
      ap1.inc(dir_);
      ket->set_member(ap1, part, pos_a);
      children_[nchildren_++] = IntType::Instance(*bra, *ket, aux, oper);
      ket->set_member(a, part, pos_a);
    }

    children_[nchildren_++] = IntType::Instance(*bra, *ket, aux, oper);

    if (!is_simple()) {  // treatment of derivative terms differs for shell sets
                         // and integrals since in computing shell sets
                         // transfer/build will occur in all 3 directions change
                         // in up to all three derivative indices will occur

      for (unsigned int xyz = 0; xyz < 3; ++xyz) {
        // is a differentiated in this direction? add another term
        if (a.deriv().d(xyz) > 0) {
          F a_der_m1(a);
          a_der_m1.deriv().dec(xyz);
          ket->set_member(a_der_m1, part, pos_a);
          children_[nchildren_++] = IntType::Instance(*bra, *ket, aux, oper);
          ket->set_member(a, part, pos_a);
        }
        // is b differentiated in this direction? add another term
        if (bm1.deriv().d(xyz) > 0) {
          F bm1_der_m1(bm1);
          bm1_der_m1.deriv().dec(xyz);
          bra->set_member(bm1_der_m1, part, pos_b);
          children_[nchildren_++] = IntType::Instance(*bra, *ket, aux, oper);
          bra->set_member(bm1, part, pos_b);
        }
      }
    }

    if (is_simple()) {
      // ( b | a ) = ( b-1_i | a+1_i )  - AB_i ( b-1_i | a ) = ( b-1_i | a+1_i )
      // + BA_i ( b-1_i | a ) the latter is amenable to generate fmadd
      // instruction
      expr_ = children_[0] + prefactors.Y_X[part][dir] * children_[1];

      // is a differentiated in this direction? add another term
      const bool aderiv = a.deriv().d(dir_) > 0;
      if (aderiv) {
        F a_der_m1(a);
        a_der_m1.deriv().dec(dir_);
        ket->set_member(a_der_m1, part, pos_a);
        children_[nchildren_++] = IntType::Instance(*bra, *ket, aux, oper);
        ket->set_member(a, part, pos_a);
      }

      // is b differentiated in this direction? add another term
      const bool bderiv = bm1.deriv().d(dir_) > 0;
      if (bderiv) {
        F bm1_der_m1(bm1);
        bm1_der_m1.deriv().dec(dir_);
        bra->set_member(bm1_der_m1, part, pos_b);
        children_[nchildren_++] = IntType::Instance(*bra, *ket, aux, oper);
        bra->set_member(bm1, part, pos_b);
      }

      if (aderiv) expr_ += Vector(a.deriv())[dir_] * children_[2];
      if (bderiv) expr_ -= Vector(b.deriv())[dir_] * children_[aderiv ? 3 : 2];
    }
  }  // a in ket, b in bra

  if (loc_a == InBra && loc_b == InKet) {
    F a(bra->member(part, pos_a));
    F b(ket->member(part, pos_b));

    // See if b-1 exists
    F bm1(b);
    bm1.dec(dir_);
    if (!exists(bm1)) {
      delete bra;
      delete ket;
      return;
    }
    ket->set_member(bm1, part, pos_b);  // set b permanently to b-1_i

    {
      F ap1(a);
      ap1.inc(dir_);
      bra->set_member(ap1, part, pos_a);
      children_[nchildren_++] = IntType::Instance(*bra, *ket, aux, oper);
      bra->set_member(a, part, pos_a);
    }

    children_[nchildren_++] = IntType::Instance(*bra, *ket, aux, oper);

    if (!is_simple()) {  // treatment of derivative terms differs for shell sets
                         // and integrals since in computing shell sets
                         // transfer/build will occur in all 3 directions change
                         // in up to all three derivative indices will occur

      for (unsigned int xyz = 0; xyz < 3; ++xyz) {
        // is a differentiated in this direction? add another term
        if (a.deriv().d(xyz) > 0) {
          F a_der_m1(a);
          a_der_m1.deriv().dec(xyz);
          bra->set_member(a_der_m1, part, pos_a);
          children_[nchildren_++] = IntType::Instance(*bra, *ket, aux, oper);
          bra->set_member(a, part, pos_a);
        }
        // is b differentiated in this direction? add another term
        if (bm1.deriv().d(xyz) > 0) {
          F bm1_der_m1(bm1);
          bm1_der_m1.deriv().dec(xyz);
          ket->set_member(bm1_der_m1, part, pos_b);
          children_[nchildren_++] = IntType::Instance(*bra, *ket, aux, oper);
          ket->set_member(bm1, part, pos_b);
        }
      }
    }

    if (is_simple()) {
      // ( a | b) = ( a+1_i | b-1_i )  + AB_i ( a | b-1_i )
      expr_ = children_[0] + prefactors.X_Y[part][dir] * children_[1];

      // is a differentiated in this direction? add another term
      const bool aderiv = a.deriv().d(dir_) > 0;
      if (aderiv) {
        F a_der_m1(a);
        a_der_m1.deriv().dec(dir_);
        bra->set_member(a_der_m1, part, pos_a);
        children_[nchildren_++] = IntType::Instance(*bra, *ket, aux, oper);
        bra->set_member(a, part, pos_a);
      }

      // is b differentiated in this direction? add another term
      const bool bderiv = bm1.deriv().d(dir_) > 0;
      if (bderiv) {
        F bm1_der_m1(bm1);
        bm1_der_m1.deriv().dec(dir_);
        ket->set_member(bm1_der_m1, part, pos_b);
        children_[nchildren_++] = IntType::Instance(*bra, *ket, aux, oper);
        ket->set_member(bm1, part, pos_b);
      }

      if (aderiv) expr_ += Vector(a.deriv())[dir_] * children_[2];
      if (bderiv) expr_ -= Vector(b.deriv())[dir_] * children_[aderiv ? 3 : 2];
    }

  }  // a in bra, b in ket

  delete bra;
  delete ket;
}

template <class IntType, class F, int part, FunctionPosition loc_a,
          unsigned int pos_a, FunctionPosition loc_b, unsigned int pos_b>
bool HRR<IntType, F, part, loc_a, pos_a, loc_b, pos_b>::register_with_rrstack()
    const {
  // This is an ugly hack -- add HRR's to the RRStack preemptively for targets
  // in which all functions not involved in transfer have zero quanta. The
  // reason is that for the HRR quartet level code to work correctly I must use
  // these particular instances of HRR to generate the source.

  using std::swap;

  // only register RRs with for shell sets
  if (TrivialBFSet<F>::result) return false;
  typedef typename IntType::BraType IBraType;
  typedef typename IntType::KetType IKetType;
  const IBraType& bra = target_->bra();
  const IKetType& ket = target_->ket();

  // check for nonzero quanta for all particles other than part
  bool nonzero_quanta = false;
  unsigned const int npart = IntType::OperatorType::Properties::np;
  for (unsigned int p = 0; p < npart; p++) {
    if (p == part) continue;
    int nfbra = bra.num_members(p);
    assert(nfbra == 1);
    for (int f = 0; f < nfbra; f++)
      if (!bra.member(p, f).zero() || !bra.member(p, f).deriv().zero())
        nonzero_quanta = true;
    int nfket = ket.num_members(p);
    assert(nfket == 1);
    for (int f = 0; f < nfket; f++)
      if (!ket.member(p, f).zero() || !ket.member(p, f).deriv().zero())
        nonzero_quanta = true;
  }
  // if all bfsets not involved in transfer have zero quanta then this instance
  // needs to be added to the stack
  if (!nonzero_quanta) {
    std::shared_ptr<RRStack> rrstack = RRStack::Instance();
    std::shared_ptr<ThisType> this_ptr =
        std::const_pointer_cast<ThisType, const ThisType>(
            std::static_pointer_cast<const ThisType, const ParentType>(
                std::enable_shared_from_this<ParentType>::shared_from_this()));
    rrstack->find(this_ptr);
    return true;
  }

  //
  // else create the needed instance of HRR
  //

  // zero out unneeded bfs'
  IBraType bra_zero(bra);
  IKetType ket_zero(ket);
  for (unsigned int p = 0; p < npart; p++) {
    if (p == part) continue;
    int nfbra = bra_zero.num_members(p);
    for (int f = 0; f < nfbra; f++) {
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
      typedef typename IKetType::bfs_type bfs_type;
      typedef typename IKetType::bfs_ref bfs_ref;
      bfs_ref bfs = ket_zero.member(p, f);
      if (!bfs.zero() || !bfs.deriv().zero()) {
        bfs_type null_bfs;
        swap(bfs, null_bfs);
      }
    }
  }

  // create a generic GenIntegralSet over a multiplicative operator
  typedef GenOper<GenMultSymmOper_Descr<IntType::OperatorType::Properties::np> >
      DummyOper;
  typedef EmptySet DummyQuanta;
  typedef GenIntegralSet<DummyOper, IncableBFSet, IBraType, IKetType,
                         DummyQuanta>
      DummyIntegral;
  DummyOper dummy_oper;
  DummyQuanta dummy_quanta(std::vector<int>(0, 0));
  std::shared_ptr<DummyIntegral> dummy_integral =
      DummyIntegral::Instance(bra_zero, ket_zero, dummy_quanta, dummy_oper);

  // Construct generic HRR and add it to the stack instead of this HRR
  typedef HRR<DummyIntegral, F, part, loc_a, pos_a, loc_b, pos_b> DummyHRR;
  std::shared_ptr<DummyHRR> dummy_hrr =
      DummyHRR::Instance(dummy_integral, dir_);
  std::shared_ptr<RRStack> rrstack = RRStack::Instance();
  rrstack->find(dummy_hrr);
  return true;
}

template <class IntType, class F, int part, FunctionPosition loc_a,
          unsigned int pos_a, FunctionPosition loc_b, unsigned int pos_b>
HRR<IntType, F, part, loc_a, pos_a, loc_b, pos_b>::~HRR() {
  oper_checks();
}

template <class IntType, class F, int part, FunctionPosition loc_a,
          unsigned int pos_a, FunctionPosition loc_b, unsigned int pos_b>
void HRR<IntType, F, part, loc_a, pos_a, loc_b, pos_b>::oper_checks() const {
  //
  // Here we check basic HRR applicability requirements on the integral class
  //

#if CHECK_SAFETY
  // part is within the range
  typedef typename IntType::OperatorType Oper;
  if (part < 0 || part >= Oper::Properties::np) {
    assert(false);
  }

  // Cannot apply when a and b are the same
  if (loc_a == loc_b && pos_a == pos_b) {
    assert(false);
  }
#endif
}

template <class IntType, class F, int part, FunctionPosition loc_a,
          unsigned int pos_a, FunctionPosition loc_b, unsigned int pos_b>
std::shared_ptr<
    typename HRR<IntType, F, part, loc_a, pos_a, loc_b, pos_b>::ChildType>
HRR<IntType, F, part, loc_a, pos_a, loc_b, pos_b>::child(unsigned int i) const {
  assert(i >= 0 && i < nchildren_);

  unsigned int nc = 0;
  for (unsigned int c = 0; c < max_nchildren_; c++) {
    if (children_[c]) {
      if (nc == i) return children_[c];
      nc++;
    }
  }
  abort();  // unreachable
}

template <class IntType, class F, int part, FunctionPosition loc_a,
          unsigned int pos_a, FunctionPosition loc_b, unsigned int pos_b>
std::string HRR<IntType, F, part, loc_a, pos_a, loc_b, pos_b>::generate_label()
    const {
  std::ostringstream os;

  os << "HRR Part " << part << " " << (loc_a == InBra ? "bra" : "ket") << " "
     << pos_a << "  " << (loc_b == InBra ? "bra" : "ket") << " " << pos_b
     << " ";

  if (loc_a == InBra) {
    F sh_a(target_->bra(part, pos_a));
    // HRR works for contracted and uncontracted non-differentiated functions
    // thus only need to create the uncontracted instances
    sh_a.uncontract();
    os << sh_a.label() << " ";

    if (loc_b == InBra) {
      F sh_b(target_->bra(part, pos_b));
      sh_b.uncontract();
      os << sh_b.label();
    } else {
      F sh_b(target_->ket(part, pos_b));
      sh_b.uncontract();
      os << sh_b.label();
    }
  } else {
    F sh_a(target_->ket(part, pos_a));
    sh_a.uncontract();
    os << sh_a.label() << " ";

    if (loc_b == InBra) {
      F sh_b(target_->bra(part, pos_b));
      sh_b.uncontract();
      os << sh_b.label();
    } else {
      F sh_b(target_->ket(part, pos_b));
      sh_b.uncontract();
      os << sh_b.label();
    }
  }

  return os.str();
}

template <class IntType, class F, int part, FunctionPosition loc_a,
          unsigned int pos_a, FunctionPosition loc_b, unsigned int pos_b>
std::string HRR<IntType, F, part, loc_a, pos_a, loc_b, pos_b>::spfunction_call(
    const std::shared_ptr<CodeContext>& context,
    const std::shared_ptr<ImplicitDimensions>& dims) const {
  std::ostringstream os;
  os << context->label_to_name(
            label_to_funcname(context->cparams()->api_prefix() + label()))
     // First argument is the library object
     << "(inteval, "
     // Second is the target
     << context->value_to_pointer(rr_target()->symbol());
  // then come children
  const unsigned int nchildren = num_children();
  for (unsigned int c = 0; c < nchildren; c++) {
    os << ", " << context->value_to_pointer(rr_child(c)->symbol());
  }
  // then dimensions of basis function sets not involved in the transfer
  unsigned int hsr = 1;
  // a cleaner way to count the number of function sets referring
  // to some particles is to construct a dummy integral and
  // use subiterator policy
  // WARNING !!!
  for (int p = 0; p < part; p++) {
    unsigned int nbra = target_->bra().num_members(p);
    for (unsigned int i = 0; i < nbra; i++) {
      SubIterator* iter = target_->bra().member_subiter(p, i);
      hsr *= iter->num_iter();
      delete iter;
    }
    unsigned int nket = target_->ket().num_members(p);
    for (unsigned int i = 0; i < nket; i++) {
      SubIterator* iter = target_->ket().member_subiter(p, i);
      hsr *= iter->num_iter();
      delete iter;
    }
  }
  // Use TaskParameters to keep track of maximum hsr
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current().params()->max_hrr_hsrank(hsr);

  // can only do a simple bra->ket or ket->bra transfer so far
  // unsigned int isr = 1;
  if (loc_a == loc_b && pos_a != 0 && pos_b != 0)
    throw CodeDoesNotExist(
        "HRR::spfunction_call -- has not been generalized yet");

  /// WARNING !!!
  unsigned int lsr = 1;
  unsigned int np = IntType::OperType::Properties::np;
  for (unsigned int p = part + 1; p < np; p++) {
    unsigned int nbra = target_->bra().num_members(p);
    for (unsigned int i = 0; i < nbra; i++) {
      SubIterator* iter = target_->bra().member_subiter(p, i);
      lsr *= iter->num_iter();
      delete iter;
    }
    unsigned int nket = target_->ket().num_members(p);
    for (unsigned int i = 0; i < nket; i++) {
      SubIterator* iter = target_->ket().member_subiter(p, i);
      lsr *= iter->num_iter();
      delete iter;
    }
  }
  // Use TaskParameters to keep track of maximum hsr
  taskmgr.current().params()->max_hrr_hsrank(hsr);

  if (expl_high_dim()) os << "," << hsr;
  if (expl_low_dim()) os << "," << lsr;
  os << ")" << context->end_of_stat() << std::endl;
  return os.str();
}

template <class IntType, class F, int part, FunctionPosition loc_a,
          unsigned int pos_a, FunctionPosition loc_b, unsigned int pos_b>
bool HRR<IntType, F, part, loc_a, pos_a, loc_b, pos_b>::expl_high_dim() const {
  bool high = true;
  if (part == 0) high = false;
  return high;
}

template <class IntType, class F, int part, FunctionPosition loc_a,
          unsigned int pos_a, FunctionPosition loc_b, unsigned int pos_b>
bool HRR<IntType, F, part, loc_a, pos_a, loc_b, pos_b>::expl_low_dim() const {
  unsigned int np = IntType::OperType::Properties::np;
  bool low = true;
  if (part == np - 1) low = false;
  // corner case: # of particles == 1
  if (np == 1) {  // to match the interface of np != 1 need to add insert dummy
                  // argument
    low = true;
  }
  return low;
}

template <class IntType, class F, int part, FunctionPosition loc_a,
          unsigned int pos_a, FunctionPosition loc_b, unsigned int pos_b>
std::shared_ptr<ImplicitDimensions>
HRR<IntType, F, part, loc_a, pos_a, loc_b, pos_b>::adapt_dims_(
    const std::shared_ptr<ImplicitDimensions>& dims) const {
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
