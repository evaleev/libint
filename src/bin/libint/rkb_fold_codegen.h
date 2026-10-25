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
#include <prefactors.h>
#include <task.h>

#include <array>
#include <cstdlib>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace libint2 {

// Shared machinery of the σ·p "element-wise fold" recurrence relations:
// CR_11_Coulombσp_11, CR_11_Coulombσpσp_11, CR_11_σpCoulombσp_11 and
// CR_11_σpσpCoulombσpσp_11. Each decomposes its target into a fixed linear
// combination of first-derivative TwoPRep children, evaluated element-wise over
// the (spectator) shell dimensions. Every such relation is defined by one table
// of terms per component (RKBFoldTerms): the constructor creates the children
// from it (rkb_fold_make_children) and generate_code() emits the element-wise
// expression from it (rkb_fold_generate_code_from_terms), so the order of the
// children and the signs in the generated code cannot disagree. The call-site
// emission, dimension adaptation and per-function scaffold are identical across
// the four relations and live here as well.

/// One term of an element-wise σ·p fold: @c coef times the TwoPRep integral
/// over copies of the functions (a, b, c, d) = (bra1, ket1, bra2, ket2) of the
/// target that carry one extra derivative quantum in Cartesian direction
/// @c deriv[0..3] (-1: no derivative on that function). The derivatives are
/// with respect to the centers of the functions (libint's convention), see
/// the rkb_operators page in engine.h for the resulting prefactors.
struct RKBFoldTerm {
  int coef;
  std::array<int, 4> deriv;
};
typedef std::vector<RKBFoldTerm> RKBFoldTerms;

/// One term of a Pauli component of the product (σ·∇)(σ·∇') of two gradients:
/// @c coef times ∂_i ∂'_j
struct RKBFoldPairTerm {
  int coef;
  int i;
  int j;
};

/// @return the terms of Pauli component @p q of
/// (σ·∇)(σ·∇') = ∇·∇' + i σ·(∇ × ∇'): q = 0 is the dot product (δ_ij),
/// q = 1..3 the x, y, z component of the cross product (the Levi-Civita symbol
/// ε_{qij}). The imaginary unit and the Pauli matrix of the cross-product
/// components are left to the user of the library (see the rkb_operators page
/// in engine.h).
inline const std::vector<RKBFoldPairTerm>& rkb_fold_pair_terms(int q) {
  static const std::vector<RKBFoldPairTerm> terms[4] = {
      {{1, 0, 0}, {1, 1, 1}, {1, 2, 2}},  // ∇·∇' = ∂_x∂'_x + ∂_y∂'_y + ∂_z∂'_z
      {{1, 1, 2}, {-1, 2, 1}},  // (∇×∇')_x = ∂_y∂'_z - ∂_z∂'_y
      {{1, 2, 0}, {-1, 0, 2}},  // (∇×∇')_y = ∂_z∂'_x - ∂_x∂'_z
      {{1, 0, 1}, {-1, 1, 0}}   // (∇×∇')_z = ∂_x∂'_y - ∂_y∂'_x
  };
  if (q < 0 || q > 3)
    throw std::invalid_argument(
        "rkb_fold_pair_terms: invalid Pauli component index " +
        std::to_string(q));
  return terms[q];
}

/// Constructor body of a fold relation: creates one TwoPRep child per term of
/// @p terms, in the order of @p terms (the order generate_code() relies on),
/// and, if @p simple, the symbolic expression Σ_t coef_t · child_t in @p expr
/// and its flop count in @p nflops. The children are created even when not
/// @p simple so that they register on the graph.
/// @tparam Factory ChildFactory<RR, GenIntegralSet_11_11<F, TwoPRep, mType>>
/// @tparam ExprPtr type of RecurrenceRelation::expr_
/// @param a,b,c,d the functions of the target: bra1, ket1, bra2, ket2
template <typename Factory, typename F, typename ExprPtr>
void rkb_fold_make_children(Factory& factory, const F& a, const F& b,
                            const F& c, const F& d, const RKBFoldTerms& terms,
                            bool simple, ExprPtr& expr, unsigned int& nflops) {
  using namespace libint2::algebra;
  using libint2::algebra::operator*;
  using libint2::prefactor::Scalar;

  const typename Factory::AuxIndexType zero_m(0u);
  std::shared_ptr<DGVertex> first_child;
  for (std::size_t t = 0; t != terms.size(); ++t) {
    const RKBFoldTerm& term = terms[t];
    F fa{a};
    F fb{b};
    F fc{c};
    F fd{d};
    if (term.deriv[0] >= 0) fa.deriv().inc(term.deriv[0]);
    if (term.deriv[1] >= 0) fb.deriv().inc(term.deriv[1]);
    if (term.deriv[2] >= 0) fc.deriv().inc(term.deriv[2]);
    if (term.deriv[3] >= 0) fd.deriv().inc(term.deriv[3]);
    const std::shared_ptr<DGVertex> child =
        factory.make_child(fa, fb, fc, fd, zero_m);
    if (!simple) continue;

    if (t == 0) {
      first_child = child;
      // a bare child is not an expression: a leading +child is combined with
      // the second term below, any other leading term is Scalar(coef) * child
      if (term.coef != 1 || terms.size() == 1) expr = Scalar(term.coef) * child;
    } else {
      const std::shared_ptr<DGVertex> lhs =
          (t == 1 && terms[0].coef == 1) ? first_child
                                         : std::shared_ptr<DGVertex>(expr);
      if (term.coef == 1)
        expr = lhs + child;
      else if (term.coef == -1)
        expr = lhs - child;
      else
        expr = lhs + Scalar(term.coef) * child;
      nflops += 1;
    }
    if (std::abs(term.coef) != 1) nflops += 1;
  }
}

/// @return the right-hand side of the element-wise fold of @p terms with child
/// t of the relation named src<t>, e.g. "src0[hsi] - src1[hsi]" or
/// "2.0*src0[hsi] - src1[hsi] - src2[hsi]"; @p nflops receives its flop count
inline std::string rkb_fold_rhs(const RKBFoldTerms& terms,
                                unsigned int& nflops) {
  std::ostringstream os;
  nflops = 0;
  for (std::size_t t = 0; t != terms.size(); ++t) {
    const int coef = terms[t].coef;
    if (t == 0) {
      if (coef < 0) os << "-";
    } else {
      os << (coef < 0 ? " - " : " + ");
      ++nflops;
    }
    if (std::abs(coef) != 1) {
      os << std::abs(coef) << ".0*";
      ++nflops;
    }
    os << "src" << t << "[hsi]";
  }
  return os.str();
}

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
  os << ")" << context->end_of_stat() << std::endl;
  return os.str();
}

/// Emits the full element-wise function: the caller supplies the per-target
/// right-hand-side expression @p rhs (e.g. "src0[hsi] + src1[hsi]") and its
/// flop count @p nflops; everything else (symbol assignment, function
/// declaration, and the element loop scaffold) is identical across the four
/// operators and emitted here. Prefer rkb_fold_generate_code_from_terms().
template <typename RR>
void rkb_fold_generate_code(const RR& self,
                            const std::shared_ptr<CodeContext>& context,
                            const std::shared_ptr<ImplicitDimensions>& dims,
                            const std::string& funcname, std::ostream& decl,
                            std::ostream& def, const std::string& rhs,
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
  def << "#include <" << context->label_to_filename(funcname) << ".h>\n\n";
  def << context->code_prefix();
  def << func_decl << context->open_block() << std::endl;
  def << context->std_function_header();
  def << "#ifdef __INTEL_COMPILER\n#pragma ivdep\n#endif\n";
  // the fold is element-wise, hence loop over all elements, including the
  // low and vector dimensions
  std::string nelem = localdims->high_label();
  if (!localdims->low_is_static() || localdims->low_label() != "1")
    nelem += "*" + localdims->low_label();
  if (!localdims->vecdim_is_static())
    nelem += "*veclen";
  else if (localdims->vecdim_label() != "1")
    nelem += "*" + localdims->vecdim_label();
  def << "for(int hsi = 0; hsi<" << nelem << "; hsi++) {\n";
  def << "target[hsi] = " << rhs << ";\n}\n";
  def << "/** Number of flops = " << nflops << " */\n";
  def << context->close_block() << std::endl;
  def << context->code_postfix();
}

/// generate_code body of a fold relation whose children were created by
/// rkb_fold_make_children() from the same @p terms: emits
/// target = Σ_t coef_t · src<t>. @p rr_name is used in the error message only.
template <typename RR>
void rkb_fold_generate_code_from_terms(
    const RR& self, const RKBFoldTerms& terms, const std::string& rr_name,
    const std::shared_ptr<CodeContext>& context,
    const std::shared_ptr<ImplicitDimensions>& dims,
    const std::string& funcname, std::ostream& decl, std::ostream& def) {
  // add_child() deduplicates children by instance (see generic_rr.h), so a
  // collision among the children would silently shift the src slots; fail the
  // export loudly instead of emitting a wrong expression
  if (self.num_children() != terms.size())
    throw std::logic_error(rr_name + "::generate_code: expected " +
                           std::to_string(terms.size()) + " children, got " +
                           std::to_string(self.num_children()));
  unsigned int nflops = 0;
  const std::string rhs = rkb_fold_rhs(terms, nflops);
  rkb_fold_generate_code(self, context, dims, funcname, decl, def, rhs, nflops);
}

}  // namespace libint2

#endif  // LIBINT_RKB_FOLD_CODEGEN_H
