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

#ifndef LIBINT_RKB_FOLD_CODEGEN_H
#define LIBINT_RKB_FOLD_CODEGEN_H

#include <dims.h>
#include <entity.h>
#include <generic_rr.h>
#include <task.h>

#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace libint2 {

// Shared code-emission helpers for the σ·p "element-wise fold" recurrence
// relations: CR_11_Coulombσp_11, CR_11_Coulombσpσp_11, CR_11_σpCoulombσp_11 and
// CR_11_σpσpCoulombσpσp_11. Each decomposes its target into a fixed linear
// combination of first-derivative TwoPRep children, evaluated element-wise over
// the (spectator) shell dimensions. Across all four, only the operator label
// and the right-hand-side expression differ; the call-site emission, dimension
// adaptation, and per-function scaffold are identical and live here so they are
// written (and corrected) once rather than copy-pasted four times.

/// adapt_dims_ body: replace the high dimension with the runtime "highdim"
/// symbol (the fold loops over the high/spectator dimension only).
inline std::shared_ptr<ImplicitDimensions> rkb_fold_adapt_dims(
    const std::shared_ptr<ImplicitDimensions>& dims) {
  auto high_dim = std::make_shared<RTimeEntity<EntityTypes::Int>>("highdim");
  return std::make_shared<ImplicitDimensions>(high_dim, dims->low(),
                                              dims->vecdim());
}

/// spfunction_call body: emit the call to the shared element-wise function,
/// passing target, the children, and total_dim (= product of all four spectator
/// shell dimensions). @p self supplies the public RR interface; @p tgt is the
/// typed target (for bra()/ket() member iteration).
template <typename RR, typename TypedTarget>
std::string rkb_fold_spfunction_call(
    const RR& self, const TypedTarget& tgt,
    const std::shared_ptr<CodeContext>& context) {
  std::ostringstream os;
  os << context->label_to_function_name(self.label()) << "(inteval, "
     << context->value_to_pointer(self.rr_target()->symbol());
  const unsigned int nc = self.num_children();
  for (unsigned int c = 0; c < nc; c++) {
    os << ", " << context->value_to_pointer(self.rr_child(c)->symbol());
  }
  // total_dim = product of all shell dims (all 4 shells are spectators)
  unsigned int total_dim = 1;
  for (unsigned int p = 0; p < 2; p++) {
    SubIterator* si = tgt->bra().member_subiter(p, 0);
    total_dim *= si->num_iter();
    delete si;
    si = tgt->ket().member_subiter(p, 0);
    total_dim *= si->num_iter();
    delete si;
  }
  os << "," << total_dim;
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  taskmgr.current().params()->max_hrr_hsrank(total_dim);
  os << ")" << context->end_of_stat() << std::endl;
  return os.str();
}

/// generate_code body: emit the full element-wise function. The caller computes
/// the per-target right-hand-side expression @p rhs (e.g.
/// "src0[hsi] + src1[hsi]") and its flop count @p nflops; everything else —
/// symbol assignment, function declaration, and the highdim loop scaffold — is
/// identical across the four operators and emitted here.
template <typename RR>
void rkb_fold_generate_code(
    const RR& self, const std::shared_ptr<CodeContext>& context,
    const std::shared_ptr<ImplicitDimensions>& dims, const std::string& funcname,
    std::ostream& decl, std::ostream& def, const std::string& rhs,
    unsigned int nflops) {
  // declare_function lives in dg.cc
  extern std::string declare_function(
      const std::shared_ptr<CodeContext>& context,
      const std::shared_ptr<ImplicitDimensions>& dims,
      const std::shared_ptr<CodeSymbols>& args, const std::string& tlabel,
      const std::string& function_descr, std::ostream& decl);

  std::shared_ptr<ImplicitDimensions> localdims = rkb_fold_adapt_dims(dims);
  // inline assign_symbols_: set symbol names on target/children and populate
  // CodeSymbols
  std::shared_ptr<CodeSymbols> symbols(new CodeSymbols);
  self.rr_target()->set_symbol("target");
  symbols->append_symbol("target");
  for (unsigned int c = 0; c < self.num_children(); c++) {
    std::string symb = "src" + std::to_string(c);
    self.rr_child(c)->set_symbol(symb);
    symbols->append_symbol(symb);
  }
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  const std::string tlabel = taskmgr.current().label();
  const std::string func_decl =
      declare_function(context, localdims, symbols, tlabel, funcname, decl);
  def << context->std_header();
  def << "#include <" << context->label_to_name(funcname) << ".h>\n\n";
  def << context->code_prefix();
  def << func_decl << context->open_block() << std::endl;
  def << context->std_function_header();
  def << "#ifdef __INTEL_COMPILER\n#pragma ivdep\n#endif\n";
  def << "for(int hsi = 0; hsi<highdim; hsi++) {\n";
  def << "target[hsi] = " << rhs << ";\n}\n";
  def << "/** Number of flops = " << nflops << " */\n";
  def << context->close_block() << std::endl;
  def << context->code_postfix();
}

}  // namespace libint2

#endif  // LIBINT_RKB_FOLD_CODEGEN_H
