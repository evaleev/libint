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

#ifndef LIBINT_COMP_11_ΣPCOULOMBΣP_11_H
#define LIBINT_COMP_11_ΣPCOULOMBΣP_11_H

#include <dims.h>
#include <entity.h>
#include <gaussoper.h>
#include <generic_rr.h>
#include <task.h>
#include <twoprep_11_11.h>

namespace libint2 {

/**
 * Computes the "Gaunt LS bilinear" integral
 *   \f$ (\mu\, \sigma\cdot\hat{p}\,\nu | 1/r_{12} | \kappa\,
 * \sigma\cdot\hat{p}\,\lambda ) \f$ in the SO(3) irreducible decomposition of
 * the rank-2 tensor \f$ T_{ab} = ( \mu \cdot \partial_a \nu | 1/r_{12} | \kappa
 * \cdot \partial_b \lambda ) \f$: 1 scalar trace + 3 antisymmetric (curl-curl)
 * + 5 symmetric-traceless = 9 components total. Each output is a small linear
 * combination of raw deriv-TwoPRep children, mirroring
 * comp_11_Coulombσpσp_11.h's pattern of trace/antisym emission.
 *
 * @tparam F basis function type. valid choices are CGShell or CGF
 */
template <typename F>
class CR_11_σpCoulombσp_11
    : public GenericRecurrenceRelation<
          CR_11_σpCoulombσp_11<F>, F,
          GenIntegralSet_11_11<F, σpCoulombσpOper, mType>> {
 public:
  typedef CR_11_σpCoulombσp_11<F> ThisType;
  typedef F BasisFunctionType;
  typedef σpCoulombσpOper OperType;
  typedef GenIntegralSet_11_11<F, σpCoulombσpOper, mType> TargetType;
  typedef GenericRecurrenceRelation<ThisType, BasisFunctionType, TargetType>
      ParentType;
  friend class GenericRecurrenceRelation<ThisType, BasisFunctionType,
                                         TargetType>;
  static const unsigned int max_nchildren = 3;

  using ParentType::Instance;

  static bool directional() { return false; }

 private:
  using ParentType::is_simple;
  using ParentType::target_;
  using ParentType::RecurrenceRelation::expr_;
  using ParentType::RecurrenceRelation::nflops_;

  /// Constructor is private, used by Instance that maintains
  /// registry of these objects
  CR_11_σpCoulombσp_11(const std::shared_ptr<TargetType>&, unsigned int = 0);

  static std::string descr() { return "CR"; }

  // --- Code sharing overrides (mirror Coulombσpσp pattern) ---
  // All shell quartets with the same quaternion component share one function.

  std::string generate_label() const override {
    return "CR_σpCoulombσp_" +
           std::to_string(target_->oper()->descr().component_index());
  }

  std::string spfunction_call(
      const std::shared_ptr<CodeContext>& context,
      const std::shared_ptr<ImplicitDimensions>& dims) const override {
    std::ostringstream os;
    os << context->label_to_function_name(this->label()) << "(inteval, "
       << context->value_to_pointer(this->rr_target()->symbol());
    const unsigned int nc = this->num_children();
    for (unsigned int c = 0; c < nc; c++) {
      os << ", " << context->value_to_pointer(this->rr_child(c)->symbol());
    }
    // total_dim = product of all shell dims (all 4 shells are spectators)
    unsigned int total_dim = 1;
    for (unsigned int p = 0; p < 2; p++) {
      SubIterator* si = target_->bra().member_subiter(p, 0);
      total_dim *= si->num_iter();
      delete si;
      si = target_->ket().member_subiter(p, 0);
      total_dim *= si->num_iter();
      delete si;
    }
    os << "," << total_dim;
    LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
    taskmgr.current().params()->max_hrr_hsrank(total_dim);
    os << ")" << context->end_of_stat() << std::endl;
    return os.str();
  }

  std::shared_ptr<ImplicitDimensions> adapt_dims_(
      const std::shared_ptr<ImplicitDimensions>& dims) const override {
    auto high_dim = std::make_shared<RTimeEntity<EntityTypes::Int>>("highdim");
    return std::make_shared<ImplicitDimensions>(high_dim, dims->low(),
                                                dims->vecdim());
  }

  /// Hand-emit the per-component irrep linear combination over deriv-ERI
  /// children. The combination depends on the target's component index.
  void generate_code(const std::shared_ptr<CodeContext>& context,
                     const std::shared_ptr<ImplicitDimensions>& dims,
                     const std::string& funcname, std::ostream& decl,
                     std::ostream& def) override {
    extern std::string declare_function(
        const std::shared_ptr<CodeContext>& context,
        const std::shared_ptr<ImplicitDimensions>& dims,
        const std::shared_ptr<CodeSymbols>& args, const std::string& tlabel,
        const std::string& function_descr, std::ostream& decl);

    std::shared_ptr<ImplicitDimensions> localdims = adapt_dims_(dims);
    std::shared_ptr<CodeSymbols> symbols(new CodeSymbols);
    this->rr_target()->set_symbol("target");
    symbols->append_symbol("target");
    for (unsigned int c = 0; c < this->num_children(); c++) {
      std::string symb = "src" + std::to_string(c);
      this->rr_child(c)->set_symbol(symb);
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

    const int comp = this->target_->oper()->descr().component_index();
    std::string rhs;
    unsigned int nflops = 0;
    switch (comp) {
      case σpCoulombσp_Descr::Scalar:
        rhs = "src0[hsi] + src1[hsi] + src2[hsi]";
        nflops = 2;
        break;
      case σpCoulombσp_Descr::AntisymX:
      case σpCoulombσp_Descr::AntisymY:
      case σpCoulombσp_Descr::AntisymZ:
      case σpCoulombσp_Descr::SymTLDiagA:
        rhs = "src0[hsi] - src1[hsi]";
        nflops = 1;
        break;
      case σpCoulombσp_Descr::SymTLDiagB:
        rhs = "2.0*src0[hsi] - src1[hsi] - src2[hsi]";
        nflops = 3;
        break;
      case σpCoulombσp_Descr::SymTLOffXY:
      case σpCoulombσp_Descr::SymTLOffXZ:
      case σpCoulombσp_Descr::SymTLOffYZ:
        rhs = "src0[hsi] + src1[hsi]";
        nflops = 1;
        break;
      default:
        throw std::runtime_error(
            "CR_11_σpCoulombσp_11::generate_code: invalid component index");
    }
    def << "target[hsi] = " << rhs << ";\n}\n";
    def << "/** Number of flops = " << nflops << " */\n";
    def << context->close_block() << std::endl;
    def << context->code_postfix();
  }
};

template <typename F>
CR_11_σpCoulombσp_11<F>::CR_11_σpCoulombσp_11(
    const std::shared_ptr<TargetType>& Tint, unsigned int)
    : ParentType(Tint, 0) {
  assert(Tint->num_func_bra(/* particle */ 0) == 1);
  assert(Tint->num_func_bra(/* particle */ 1) == 1);
  assert(Tint->num_func_ket(/* particle */ 0) == 1);
  assert(Tint->num_func_ket(/* particle */ 1) == 1);

  F a(Tint->bra(0, 0));
  F b(Tint->ket(0, 0));
  F c(Tint->bra(1, 0));
  F d(Tint->ket(1, 0));

  const auto& oper = Tint->oper();

  if (a.contracted() || b.contracted() || c.contracted() || d.contracted())
    return;

  using namespace libint2::algebra;
  using namespace libint2::prefactor;
  using libint2::algebra::operator*;

  const mType zero_m(0u);

  ChildFactory<ThisType,
               GenIntegralSet_11_11<BasisFunctionType, TwoPRep, mType>>
      factory(this);

  // Chemist notation: (a b | op c op d). σ·p acts on one function per electron
  // — direction `i` on ket(0,0) = b (electron 1), direction `j` on ket(1,0) = d
  // (electron 2). Each output is an SO(3) irrep combination of the raw 3×3
  // T_{ij} dyadic; case bodies build only the children each combination needs.
  constexpr auto x = 0;
  constexpr auto y = 1;
  constexpr auto z = 2;

  auto T = [&](int i, int j) {
    F b_d{b};
    b_d.deriv().inc(i);
    F d_d{d};
    d_d.deriv().inc(j);
    return factory.make_child(a, b_d, c, d_d, zero_m);
  };

  switch (oper->descr().component_index()) {
    case σpCoulombσp_Descr::Scalar: {
      auto Txx = T(x, x);
      auto Tyy = T(y, y);
      auto Tzz = T(z, z);
      if (is_simple()) {
        expr_ = Txx + Tyy + Tzz;
        nflops_ += 2;
      }
    } break;
    case σpCoulombσp_Descr::AntisymX: {
      auto Tyz = T(y, z);
      auto Tzy = T(z, y);
      if (is_simple()) {
        expr_ = Tyz - Tzy;
        nflops_ += 1;
      }
    } break;
    case σpCoulombσp_Descr::AntisymY: {
      auto Tzx = T(z, x);
      auto Txz = T(x, z);
      if (is_simple()) {
        expr_ = Tzx - Txz;
        nflops_ += 1;
      }
    } break;
    case σpCoulombσp_Descr::AntisymZ: {
      auto Txy = T(x, y);
      auto Tyx = T(y, x);
      if (is_simple()) {
        expr_ = Txy - Tyx;
        nflops_ += 1;
      }
    } break;
    case σpCoulombσp_Descr::SymTLDiagA: {
      auto Txx = T(x, x);
      auto Tyy = T(y, y);
      if (is_simple()) {
        expr_ = Txx - Tyy;
        nflops_ += 1;
      }
    } break;
    case σpCoulombσp_Descr::SymTLDiagB: {
      // 2·T_zz − T_xx − T_yy: child order (Tzz, Txx, Tyy) matches generate_code
      auto Tzz = T(z, z);
      auto Txx = T(x, x);
      auto Tyy = T(y, y);
      if (is_simple()) {
        expr_ = Scalar(2) * Tzz - Txx - Tyy;
        nflops_ += 3;
      }
    } break;
    case σpCoulombσp_Descr::SymTLOffXY: {
      auto Txy = T(x, y);
      auto Tyx = T(y, x);
      if (is_simple()) {
        expr_ = Txy + Tyx;
        nflops_ += 1;
      }
    } break;
    case σpCoulombσp_Descr::SymTLOffXZ: {
      auto Txz = T(x, z);
      auto Tzx = T(z, x);
      if (is_simple()) {
        expr_ = Txz + Tzx;
        nflops_ += 1;
      }
    } break;
    case σpCoulombσp_Descr::SymTLOffYZ: {
      auto Tyz = T(y, z);
      auto Tzy = T(z, y);
      if (is_simple()) {
        expr_ = Tyz + Tzy;
        nflops_ += 1;
      }
    } break;
    default:
      throw std::runtime_error("CR_11_σpCoulombσp_11: invalid component index");
  }

}  // CR_11_σpCoulombσp_11<F>::CR_11_σpCoulombσp_11

}  // namespace libint2

#endif  // LIBINT_COMP_11_ΣPCOULOMBΣP_11_H
