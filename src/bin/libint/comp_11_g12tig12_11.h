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

#ifndef _libint2_src_bin_libint_cr11g12tig1211_h_
#define _libint2_src_bin_libint_cr11g12tig1211_h_

#include <gaussoper.h>
#include <generic_rr.h>

#include "task.h"

namespace libint2 {

/** Compute relation for 2-e integrals of the G12_Ti_G12 operators.
 */
template <class BFSet>
class CR_11_G12TiG12_11 : public GenericRecurrenceRelation<
                              CR_11_G12TiG12_11<BFSet>, BFSet,
                              GenIntegralSet_11_11<BFSet, G12TiG12, mType> > {
 public:
  typedef CR_11_G12TiG12_11 ThisType;
  typedef BFSet BasisFunctionType;
  typedef GenIntegralSet_11_11<BFSet, G12TiG12, mType> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BFSet, TargetType> ParentType;
  friend class GenericRecurrenceRelation<ThisType, BFSet, TargetType>;
  static const unsigned int max_nchildren = 1;

  using ParentType::Instance;

  /// This relation is not directional
  static bool directional() { return false; }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;

  /// Constructor is private, used by ParentType::Instance that maintains
  /// registry of these objects
  CR_11_G12TiG12_11(const std::shared_ptr<TargetType>&, unsigned int dir);
  static std::string descr() { return "CR"; }

#if LIBINT_ENABLE_GENERIC_CODE
  /// Implementation of RecurrenceRelation::has_generic()
  bool has_generic(
      const std::shared_ptr<CompilationParameters>& cparams) const override;
  /// Implementation of RecurrenceRelation::generic_header()
  std::string generic_header() const override { return "GenericScale.h"; }
  /// Implementation of RecurrenceRelation::generic_instance()
  std::string generic_instance(
      const std::shared_ptr<CodeContext>& context,
      const std::shared_ptr<CodeSymbols>& args) const override;
#endif
};

template <class F>
CR_11_G12TiG12_11<F>::CR_11_G12TiG12_11(const std::shared_ptr<TargetType>& Tint,
                                        unsigned int dir)
    : ParentType(Tint, dir) {
  if (dir != 0) return;
  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using namespace libint2::braket;
  // kinetic energy of which electron?
  // const int i = target_->oper()->descr().K();
  const R12kG12 G2(2);

  F a(Tint->bra(0, 0));
  F b(Tint->ket(0, 0));
  F c(Tint->bra(1, 0));
  F d(Tint->ket(1, 0));

  if (target_->oper()->descr().contracted()) return;

  // [G12,[T1,G12]] = [G12,[T2,G12]]
  {
    typedef GenIntegralSet_11_11<BasisFunctionType, R12kG12, mType> ChildType;
    ChildFactory<ThisType, ChildType> factory(this);
    auto ab_G2_cd = factory.make_child(a, b, c, d, 0u, G2);
    if (is_simple()) {
      expr_ = Scalar("R12_2_G12_scale_to_G12T1G12") * ab_G2_cd;
      nflops_ += 1;
    }
  }
}

#if LIBINT_ENABLE_GENERIC_CODE
template <class F>
bool CR_11_G12TiG12_11<F>::has_generic(
    const std::shared_ptr<CompilationParameters>& cparams) const {
  F sh_a(target_->bra(0, 0));
  F sh_b(target_->ket(0, 0));
  F sh_c(target_->bra(1, 0));
  F sh_d(target_->ket(1, 0));
  const unsigned int max_opt_am = cparams->max_am_opt();
  // to generate optimized code for xxxx integral need to generate specialized
  // code for up to (x+x)0(x+x)0 integrals
  if (!TrivialBFSet<F>::result &&
      (sh_a.norm() > max_opt_am || sh_b.norm() > max_opt_am ||
       sh_c.norm() > max_opt_am || sh_d.norm() > max_opt_am))
    return true;
  return false;
}

template <class F>
std::string CR_11_G12TiG12_11<F>::generic_instance(
    const std::shared_ptr<CodeContext>& context,
    const std::shared_ptr<CodeSymbols>& args) const {
  std::ostringstream oss;

  const bool vec = (context->cparams()->max_vector_length() != 1);
  if (vec)
    oss << "_libint2_static_api_scale_vec_short_(";
  else
    oss << "_libint2_static_api_scale_short_(";

  const unsigned int nargs = args->n();
  for (unsigned int a = 0; a < nargs; a++) {
    oss << args->symbol(a) << ",";
  }

  oss << target_->size() << "*" << (vec ? "inteval->veclen" : "1")
      << ",inteval->R12_2_G12_scale_to_G12T1G12"
      << (vec ? ",inteval->veclen" : "[0]") << ");";

  // HACK alert ... unfortunately it's not completely possible to figure out all
  // "precomputed" symbols from the DAG alone force R12_2_G12_scale_to_G12T1G12
  // on the list of symbols
  {
    LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
    std::list<std::string> forced_symbol;
    forced_symbol.push_back(std::string("R12_2_G12_scale_to_G12T1G12"));
    taskmgr.current().symbols()->add(forced_symbol);
  }

  return oss.str();
}
#endif

};  // namespace libint2

#endif
