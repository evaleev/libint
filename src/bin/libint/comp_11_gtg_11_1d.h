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

#ifndef _libint2_src_bin_libint_cr11gtg111d_h_
#define _libint2_src_bin_libint_cr11gtg111d_h_

#include <generic_rr.h>
#include <util_types.h>

namespace libint2 {

/** Compute relation for 1-dimensional Gaussian-type geminal integrals
 */
template <CartesianAxis Axis>
class CR_11_GTG_11_1d
    : public GenericRecurrenceRelation<
          CR_11_GTG_11_1d<Axis>, CGShell1d<Axis>,
          GenIntegralSet_11_11<CGShell1d<Axis>, GTG_1d, EmptySet> > {
 public:
  typedef CR_11_GTG_11_1d ThisType;
  typedef CGShell1d<Axis> BasisFunctionType;
  typedef GenIntegralSet_11_11<BasisFunctionType, GTG_1d, EmptySet> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BasisFunctionType, TargetType>
      ParentType;
  friend class GenericRecurrenceRelation<ThisType, BasisFunctionType,
                                         TargetType>;
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
  CR_11_GTG_11_1d(const std::shared_ptr<TargetType>&, unsigned int dir);
  static std::string descr() { return "CR"; }

#if LIBINT_ENABLE_GENERIC_CODE
  /// Implementation of RecurrenceRelation::has_generic()
  bool has_generic(
      const std::shared_ptr<CompilationParameters>& cparams) const override {
    return true;
  }
  /// Implementation of RecurrenceRelation::generic_header()
  std::string generic_header() const override { return "VRR_GTG_1d_xx_xx.h"; }
  /// Implementation of RecurrenceRelation::generic_instance()
  std::string generic_instance(
      const std::shared_ptr<CodeContext>& context,
      const std::shared_ptr<CodeSymbols>& args) const override;
#endif
};

template <CartesianAxis Axis>
CR_11_GTG_11_1d<Axis>::CR_11_GTG_11_1d(const std::shared_ptr<TargetType>& Tint,
                                       unsigned int dir)
    : ParentType(Tint, dir) {
  if (dir != 0) return;
  const GTG_1d oper;

  BasisFunctionType _0(0);

  if (target_->oper()->descr().contracted()) return;

  ChildFactory<ThisType, TargetType> factory(this);
  auto _00_GTG_00 = factory.make_child(_0, _0, _0, _0, 0u, oper);
}

#if LIBINT_ENABLE_GENERIC_CODE
template <CartesianAxis Axis>
std::string CR_11_GTG_11_1d<Axis>::generic_instance(
    const std::shared_ptr<CodeContext>& context,
    const std::shared_ptr<CodeSymbols>& args) const {
  std::ostringstream oss;
  auto a = target_->bra(0, 0);
  auto b = target_->ket(0, 0);
  auto c = target_->bra(1, 0);
  auto d = target_->kt(1, 0);

  const bool vec = (context->cparams()->max_vector_length() != 1);

  oss << "VRR_GTG_1d_xx_xx::compute<" << to_string<Axis> << "," << a[0] << ","
      << b[0] << "," << c[0] << "," << d[0] << "," << (vec ? "true" : "false")
      << ">::compute(";

  const unsigned int nargs = args->n();
  for (unsigned int a = 0; a < nargs; a++) {
    oss << args->symbol(a) << ",";
  }

  oss << ",inteval->_00_GTG1d_00_" << to_string(Axis) << ");";

  // force some quantities on the list of task symbols
  {
    LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
    std::list<std::string> forced_symbol;
    forced_symbol.push_back(std::string("AB_") + to_string(Axis));
    forced_symbol.push_back(std::string("CD_") + to_string(Axis));
    forced_symbol.push_back(std::string("R12kG12_pfac0_0_") + to_string(Axis));
    forced_symbol.push_back(std::string("R12kG12_pfac0_1_") + to_string(Axis));
    forced_symbol.push_back(std::string("R12kG12_pfac1_0"));
    forced_symbol.push_back(std::string("R12kG12_pfac1_1"));
    forced_symbol.push_back(std::string("R12kG12_pfac2"));
    taskmgr.current().symbols()->add(forced_symbol);
  }

  return oss.str();
}
#endif

};  // namespace libint2

#endif
