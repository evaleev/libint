/*
 *  Copyright (C) 2004-2017 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef INCLUDE_LIBINT2_LCAO_MOLDEN_H_
#define INCLUDE_LIBINT2_LCAO_MOLDEN_H_

#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

#include <libint2/atom.h>
#include <libint2/basis.h>
#include <libint2/cgshell_ordering.h>
#include <libint2/chemistry/elements.h>
#include <libint2/shell.h>
#include <libint2/shgshell_ordering.h>

#pragma GCC diagnostic push
#pragma GCC system_header
#include <Eigen/Core>
#pragma GCC diagnostic pop

namespace libint2 {
namespace molden {

/// Exports LCAO coefficients in [Molden
/// format](http://www.cmbi.ru.nl/molden/molden_format.html).
class Export {
 public:
  /// @tparam Coeffs the type of LCAO coefficient matrix
  /// @tparam Energies the type of LCAO energy vector
  /// @tparam Occs the type of LCAO occupancy vector
  /// @param atoms the set of atoms
  /// @param basis the set of shells; must meet Molden requirements (see below)
  /// @param coefficients the matrix of LCAO coefficients (columns are LCAOs,
  ///        rows are AOs; AOs are ordered according to the order of shells in
  ///        \c basis and by the ordering conventions of this Libint
  ///        configuration)
  /// @param occupancies the vector of occupancies (size = # LCAOs)
  /// @param energies the vector of energies (size = # of LCAOs); the default is
  ///        to assign zero to each LCAO
  /// @param symmetry_labels the vector of symmetry labels (size = # LCAOs); the
  ///        default is to assign empty label to each LCAO
  /// @param spincases the vector of spin cases (size = # LCAOs; true = spin-up
  ///        or m_s=1/2, false = spin-down or m_s=-1/2); the default is
  ///        to assign spin-up to each LCAO
  /// @throw std::runtime_error if the basis does not conforms Molden
  ///        requirements
  /// @note Molden can only handle basis sets that:
  /// - p (l=1) shells are Cartesian, not solid harmonics
  /// - d, f, and g (l=2..4) shells are all Cartesian or all solid harmonics
  /// - there are no shells with l>5
  template <typename Coeffs, typename Occs, typename Energies = Eigen::MatrixXd>
  Export(const std::vector<Atom>& atoms, const std::vector<Shell>& basis,
         const Coeffs& coefficients, const Occs& occupancies,
         const Energies& energies = Energies(),
         const std::vector<std::string>& symmetry_labels =
             std::vector<std::string>(),
         const std::vector<bool>& spincases = std::vector<bool>())
      : atoms_(atoms),
        basis_(validate(basis)),
        coefs_(coefficients),
        occupancies_(occupancies),
        energies_(energies),
        labels_(symmetry_labels),
        spins_(spincases) {
    initialize_bf_map();
  }

  /// writes "[Molden Format]" to ostream \c os
  void write_prologue(std::ostream& os) const {
    os << "[Molden Format]" << std::endl;
  }

  /// writes the "[Atoms]" section to ostream \c os
  void write_atoms(std::ostream& os) const {
    os << "[Atoms] AU" << std::endl;

    os.fill(' ');
    os << std::fixed << std::setprecision(8);
    auto iatom = 0;
    for (const auto& atom : atoms_) {
      auto Z = atom.atomic_number;
      os << std::setw(4) << libint2::chemistry::element_info[Z - 1].symbol
         << std::setw(6) << (iatom + 1) << std::setw(6) << Z << std::setw(14)
         << atom.x << std::setw(14) << atom.y << std::setw(14) << atom.z
         << std::endl;
      ++iatom;
    }
  }

  /// writes the "[GTO]" section, as well as optional Cartesian/solid harmonics
  /// keywords, to ostream \c os
  void write_basis(std::ostream& os) const {
    os << "[GTO]" << std::endl;
    for (size_t iatom = 0; iatom < atoms_.size(); ++iatom) {
      os << std::setw(4) << (iatom + 1) << std::setw(4) << 0 << std::endl;
      for (auto ish : atom2shell_[iatom]) {
        const Shell& sh = basis_.at(ish);
        if (sh.contr.size() == 1) {
          const auto& contr = sh.contr[0];
          const auto nprim = contr.coeff.size();
          os << std::setw(4) << Shell::am_symbol(contr.l) << std::setw(6)
             << nprim << std::setw(6) << "1.00" << std::endl;
          for (int iprim = 0; iprim < nprim; ++iprim) {
            os << std::scientific << std::uppercase << std::setprecision(10);
            os << std::setw(20) << sh.alpha[iprim] << std::setw(20)
               << sh.coeff_normalized(0, iprim) << std::endl;
          }  // end loop over primitives
        }    // end if ncontraction == 1
        else {
          assert(false);  // Not implemented
        }
      }  // end loop over shells on center
      // format calls for a blank line here
      os << std::endl;
    }  // end loop over centers

    // write solid harmonic/cartesian tags
    {
      // default is cartesians throughout
      if (dfg_is_cart_[0]) {   // cartesian d
        if (!dfg_is_cart_[1])  // solid harmonic f
          os << "[7F]" << std::endl;
      } else {                   // solid harmonic d
        if (!dfg_is_cart_[1]) {  // solid harmonic f
          os << "[5D7F]" << std::endl;
        } else {  // cartesian f
          os << "[5D10F]" << std::endl;
        }
      }
      if (!dfg_is_cart_[2])  // solid harmonic g
        os << "[9G]" << std::endl;
    }
  }

  /// writes the "[MO]" section to ostream \c os
  void write_lcao(std::ostream& os) const {
    // TODO Open shell cases
    os << "[MO]" << std::endl;
    os << std::fixed << std::setprecision(10);
    for (int imo = 0; imo < coefs_.cols(); ++imo) {
      os << std::setw(8) << "Sym=" << (labels_.empty() ? "" : labels_.at(imo))
         << std::endl
         << std::setw(8) << "Ene=" << std::setw(16)
         << (energies_.rows() == 0 ? 0.0 : energies_(imo)) << std::endl
         << std::setw(8) << "Spin="
         << (spins_.empty() ? "Alpha" : (spins_.at(imo) ? "Alpha" : "Beta"))
         << std::endl
         << std::setw(8) << "Occup=" << occupancies_(imo) << std::endl;
      for (int iao = 0; iao < coefs_.rows(); ++iao) {
        const auto C_ao_mo = coefs_(iao, imo);
        if (std::abs(C_ao_mo) >= 5e-11) {
          os << std::setw(6) << (ao_map_[iao] + 1) << std::setw(16)
             << C_ao_mo << std::endl;
        }
      }  // end loop over AOs
    }    // end loop over MOs
  }

  /// writes "prologue", atoms, basis, and LCAOs to ostream \c os
  void write(std::ostream& os) const {
    write_prologue(os);
    write_atoms(os);
    write_basis(os);
    write_lcao(os);
  }

  /// same as write(ostream), but creates new file named \c filename
  void write(const std::string& filename) const {
    std::ofstream os(filename);
    write_prologue(os);
    write_atoms(os);
    write_basis(os);
    write_lcao(os);
  }

 private:
  const std::vector<Atom>& atoms_;
  const std::vector<Shell>& basis_;
  Eigen::MatrixXd coefs_;
  Eigen::VectorXd occupancies_;
  Eigen::VectorXd energies_;
  std::vector<std::string> labels_;
  std::vector<bool> spins_;
  mutable bool
      dfg_is_cart_[3];  // whether {d, f, g} shells are cartesian (true) or
                        // solid harmonics (false)
  std::vector<std::vector<long>>
      atom2shell_;  // maps atom -> shell indices in basis_
  std::vector<long>
      ao_map_;  // maps from AO order of basis_ to the order assumed by Molden
                // (atoms->shells, bf in shells ordered in the Molden order)

  /// @throw std::runtime_error if the basis does not conforms Molden
  ///        requirements
  const std::vector<Shell>& validate(const std::vector<Shell>& shells) const {
    bool dfg_found[] = {false, false, false};
    for (const auto& shell : shells) {
      for (const auto& contr : shell.contr) {
        if (contr.l > 4)
          throw std::runtime_error(
              "molden::Export cannot handle shells with l > 4");

        switch (contr.l) {
          case 1:
            if (contr.pure)
              throw std::runtime_error(
                  "molden::Export cannot handle solid harmonics p shells");
            break;
          case 2:
          case 3:
          case 4: {
            if (!dfg_found[contr.l - 2]) {
              dfg_is_cart_[contr.l - 2] = !contr.pure;
            }
            if (!contr.pure ^ dfg_is_cart_[contr.l - 2])
              throw std::runtime_error(
                  "molden::Export only supports all-Cartesian or "
                  "all-solid-harmonics d/f/g shells");
          }

          default: {}  // l = 0 is fine
        }
      }
    }
    return shells;
  }

  void initialize_bf_map() {
    atom2shell_ = BasisSet::atom2shell(atoms_, basis_);

    const auto nao = BasisSet::nbf(basis_);
    ao_map_.resize(nao);
    assert(nao == coefs_.rows());
    const auto shell2ao = BasisSet::compute_shell2bf(basis_);
    long ao_molden = 0;
    for (size_t iatom = 0; iatom < atoms_.size(); ++iatom) {
      for (auto ish : atom2shell_[iatom]) {
        auto ao = shell2ao[ish];  // refers to order assumed by coefs
        const auto& shell = basis_[ish];
        const auto ncontr = shell.contr.size();
        for (int c = 0; c != ncontr; ++c) {
          const auto l = shell.contr[c].l;
          const auto pure = shell.contr[c].pure;
          if (pure) {
            int m;
            FOR_SOLIDHARM_MOLDEN(l, m)
            const auto ao_in_shell = INT_SOLIDHARMINDEX(l, m);
            ao_map_[ao + ao_in_shell] = ao_molden;
            ++ao_molden;
            END_FOR_SOLIDHARM_MOLDEN
            ao += 2 * l + 1;
          } else {
            int i, j, k;
            FOR_CART_MOLDEN(i, j, k, l)
            const auto ao_in_shell = INT_CARTINDEX(l, i, j);
            ao_map_[ao + ao_in_shell] = ao_molden;
            ++ao_molden;
            END_FOR_CART_MOLDEN
            ao += INT_NCART(l);
          }
        }  // contraction loop
      }
    }
  }

};  // Export

}  // namespace molden
}  // namespace libint2

#endif  // INCLUDE_LIBINT2_LCAO_MOLDEN_H_
