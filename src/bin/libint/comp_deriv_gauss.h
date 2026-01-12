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

#ifndef _libint2_src_bin_libint_compderivgauss_h_
#define _libint2_src_bin_libint_compderivgauss_h_

#include <generic_rr.h>

#include <set>

namespace libint2 {

class CR_DerivGauss_GenericInstantiator {
  static CR_DerivGauss_GenericInstantiator instance_;

  CR_DerivGauss_GenericInstantiator();  // this is a singleton
  ~CR_DerivGauss_GenericInstantiator();

  // pairs of L,vectorize specify the instances of GenericGaussDeriv template to
  // be created
  std::set<std::pair<unsigned int, bool> > template_instances_;

 public:
  static CR_DerivGauss_GenericInstantiator& instance();
  void add(unsigned int L, bool vectorize);
};

/** Compute relation for (geometric) derivative Gaussian ints of generic type \c
 * IntType . It either reduces derivative Gaussian to other Gaussians, or
 * optionally relates derivative integral to other integrals using translational
 * invariance relation. The advantage of the latter is that it can be applied to
 * contracted and solid harmonics integrals, and reduces the total number of
 * derivatives to be computed.
 *
 * @tparam IntType integral type
 * @tparam part particle index of the function to be differentiated
 * @tparam where position of the function to be differentiated
 * @tparam trans_inv_part if non-negative, specifies the particle index of the
 * function whose derivatives will be computed using translational invariance
 * @tparam trans_inv_where the position of the function
 *         whose derivatives will be computed using translational invariance
 * @note Translational will be used iff
 *         @code part == trans_inv_part && where == trans_inv_where @endcode
 */
template <class IntType, int part, FunctionPosition where,
          int trans_inv_part = -1, FunctionPosition trans_inv_where = InBra>
class CR_DerivGauss
    : public GenericRecurrenceRelation<
          CR_DerivGauss<IntType, part, where, trans_inv_part, trans_inv_where>,
          typename IntType::BasisFunctionType, IntType> {
 private:
  static constexpr auto trans_inv_oper =
      not IntType::OperType::Properties::odep;
  // can use translational invariance iff the operator is
  // translationally-invariant
  static constexpr auto using_trans_inv =
      trans_inv_oper && (part == trans_inv_part) && (where == trans_inv_where);

 public:
  typedef CR_DerivGauss ThisType;
  typedef typename IntType::BasisFunctionType BasisFunctionType;
  typedef IntType TargetType;
  typedef GenericRecurrenceRelation<ThisType, BasisFunctionType, TargetType>
      ParentType;
  friend class GenericRecurrenceRelation<ThisType, BasisFunctionType,
                                         TargetType>;
  // # of children varies:
  // 1. translational invariance makes num_functions - 1 children
  // 2. direct differentiation always makes at most 2 Gaussians
  static constexpr auto max_nchildren =
      using_trans_inv ? (IntType::num_bf - 1) : 2u;

  using ParentType::Instance;

  /// always directional! Cartesian derivatives are applied in a particular
  /// direction
  static bool directional() { return true; }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;

  /// Constructor is private, used by ParentType::Instance that maintains
  /// registry of these objects
  CR_DerivGauss(const std::shared_ptr<TargetType>&, unsigned int dir);

  static std::string descr() { return "CR_DerivGauss"; }
  /** Re-Implementation of GenericRecurrenceRelation::generate_label():
      CR Deriv relations are independent of m (it never appears anywhere in
     equations), hence to avoid generating identical code make sure that the
     (unique) label has m=0. */
  std::string generate_label() const override {
    typedef typename TargetType::AuxIndexType mType;
    static std::shared_ptr<mType> aux0(new mType(0u));
    std::ostringstream os;
    os << descr() << "P" << part << to_string(where)
       << genintegralset_label(target_->bra(), target_->ket(), aux0,
                               target_->oper());
    return os.str();
  }

#if LIBINT_ENABLE_GENERIC_CODE
  /// Implementation of RecurrenceRelation::has_generic()
  bool has_generic(
      const std::shared_ptr<CompilationParameters>& cparams) const override;
  /// Implementation of RecurrenceRelation::generic_header()
  std::string generic_header() const override;
  /// Implementation of RecurrenceRelation::generic_instance()
  std::string generic_instance(
      const std::shared_ptr<CodeContext>& context,
      const std::shared_ptr<CodeSymbols>& args) const override;
#endif
};

template <class IntType, int part, FunctionPosition where, int trans_inv_part,
          FunctionPosition trans_inv_where>
CR_DerivGauss<IntType, part, where, trans_inv_part, trans_inv_where>::
    CR_DerivGauss(const std::shared_ptr<TargetType>& Tint, unsigned int dir)
    : ParentType(Tint, dir) {
  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using namespace libint2::braket;
  typedef BasisFunctionType F;
  const F& _1 = unit<F>(
      is_simple() ? dir : 0);  // for shell sets increment the first index

  const typename IntType::AuxQuantaType& aux = Tint->aux();
  const typename IntType::OperType& oper = Tint->oper();

  // WARNING assuming one function per position

  // the Gaussian must be differentiated in direction dir
  {
    if (where == InBra && Tint->bra(part, 0).deriv().d(dir) == 0) return;
    if (where == InKet && Tint->ket(part, 0).deriv().d(dir) == 0) return;
  }

  // if not using translational invariance ...
  // can expand derivatives of primitive Gaussians only
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
      auto int_ap1 = this->add_child(IntType::Instance(*bra, *ket, aux, oper));
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
        auto int_am1 =
            this->add_child(IntType::Instance(*bra, *ket, aux, oper));
        bra->set_member(a, part, 0);
        if (is_simple()) {
          expr_ -= Scalar(a[dir]) * int_am1;
          nflops_ += 2;
        }
      }
      return;
    }

    if (where == InKet) {
      F a(ket->member(part, 0));

      // add a+1
      F ap1(ket->member(part, 0) + _1);
      ap1.deriv().dec(dir);
      ket->set_member(ap1, part, 0);
      auto int_ap1 = this->add_child(IntType::Instance(*bra, *ket, aux, oper));
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
        auto int_am1 =
            this->add_child(IntType::Instance(*bra, *ket, aux, oper));
        ket->set_member(a, part, 0);
        if (is_simple()) {
          expr_ -= Scalar(a[dir]) * int_am1;
          nflops_ += 2;
        }
      }
      return;
    }
  } else {  // use translational invariance

    // remove one deriv quantum from the target function
    if (where == InBra) bra->member(part, 0).deriv().dec(dir);
    if (where == InKet) ket->member(part, 0).deriv().dec(dir);

    int term_count = 0;
    for (int p = 0; p != IntType::num_particles; ++p) {
      if (p != trans_inv_part || trans_inv_where != InBra) {
        F a(bra->member(p, 0));
        if (not a.is_unit()) {
          F da(a);
          da.deriv().inc(dir);
          bra->set_member(da, p, 0);
          auto int_da =
              this->add_child(IntType::Instance(*bra, *ket, aux, oper));
          bra->set_member(a, p, 0);
          if (is_simple()) {
            std::ostringstream oss;
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
          auto int_da =
              this->add_child(IntType::Instance(*bra, *ket, aux, oper));
          ket->set_member(a, p, 0);
          if (is_simple()) {
            std::ostringstream oss;
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

  return;
}

#if LIBINT_ENABLE_GENERIC_CODE
template <class IntType, int part, FunctionPosition where, int trans_inv_part,
          FunctionPosition trans_inv_where>
bool CR_DerivGauss<IntType, part, where, trans_inv_part, trans_inv_where>::
    has_generic(const std::shared_ptr<CompilationParameters>& cparams) const {
  // not implemented generic relation for translational invariance yet
  if (using_trans_inv) return false;

  if (TrivialBFSet<BasisFunctionType>::result) return false;

  // generate generic code if the average quantum number > max_l_opt
  const unsigned int max_opt_am = cparams->max_am_opt();
  unsigned int am_total = 0;
  unsigned int nfunctions = 0;
  const unsigned int np = IntType::OperType::Properties::np;
  for (unsigned int p = 0; p < np; p++) {
    unsigned int nbra = target_->bra().num_members(p);
    for (unsigned int i = 0; i < nbra; i++) {
      am_total += target_->bra(p, i).norm();
      ++nfunctions;
    }
    unsigned int nket = target_->ket().num_members(p);
    for (unsigned int i = 0; i < nket; i++) {
      am_total += target_->ket(p, i).norm();
      ++nfunctions;
    }
  }
  if (am_total > max_opt_am * nfunctions) return true;

  // else generate explicit code
  return false;
}

template <class IntType, int part, FunctionPosition where, int trans_inv_part,
          FunctionPosition trans_inv_where>
std::string CR_DerivGauss<IntType, part, where, trans_inv_part,
                          trans_inv_where>::generic_header() const {
  return std::string("GenericGaussDeriv.h");
}

template <class IntType, int part, FunctionPosition where, int trans_inv_part,
          FunctionPosition trans_inv_where>
std::string
CR_DerivGauss<IntType, part, where, trans_inv_part, trans_inv_where>::
    generic_instance(const std::shared_ptr<CodeContext>& context,
                     const std::shared_ptr<CodeSymbols>& args) const {
  std::ostringstream oss;

  oss << "using namespace libint2;" << std::endl;

  BasisFunctionType sh(where == InBra ? target_->bra(part, 0)
                                      : target_->ket(part, 0));

  const unsigned int L = sh.norm();
  const bool vectorize =
      (context->cparams()->max_vector_length() == 1) ? false : true;
  oss << "libint2::GenericGaussDeriv<" << L << ","
      << (vectorize ? "true" : "false") << ">::compute(inteval";

  oss << "," << args->symbol(0);  // target
  const unsigned int nargs = args->n();
  for (unsigned int a = 1; a < nargs; a++) {
    oss << "," << args->symbol(a);
  }
  // L == 0 => second argument not needed
  if (nargs == 2) oss << ",0";

  // then dimensions of basis function sets not involved in the transfer
  unsigned int hsr = 1;
  unsigned int lsr = 1;
  const unsigned int np = IntType::OperType::Properties::np;
  // a cleaner way to count the number of function sets referring
  // to some particles is to construct a dummy integral and
  // use subiterator policy
  // WARNING !!!
  for (int p = 0; p < static_cast<int>(np); p++) {
    unsigned int nbra = target_->bra().num_members(p);
    assert(nbra == 1);
    for (unsigned int i = 0; i < nbra; i++) {
      SubIterator* iter = target_->bra().member_subiter(p, i);
      if (p < part || (p == part && where == InKet)) hsr *= iter->num_iter();
      // skip p == part && where == InBra
      if (p > part) lsr *= iter->num_iter();
      delete iter;
    }
    unsigned int nket = target_->ket().num_members(p);
    assert(nket == 1);
    for (unsigned int i = 0; i < nket; i++) {
      SubIterator* iter = target_->ket().member_subiter(p, i);
      if (p < part) hsr *= iter->num_iter();
      // skip p == part && where == InKet
      if (p > part || (p == part && where == InBra)) lsr *= iter->num_iter();
      delete iter;
    }
  }
  oss << "," << hsr << "," << lsr;

  // direction
  oss << "," << this->dir();

  // two_alpha prefactor
  oss << ",inteval->two_alpha" << part << "_"
      << (where == InBra ? "bra" : "ket");

  oss << ");";

  CR_DerivGauss_GenericInstantiator::instance().add(L, vectorize);

  return oss.str();
}
#endif  // #if !LIBINT_ENABLE_GENERIC_CODE

};  // namespace libint2

#endif
