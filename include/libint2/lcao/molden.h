/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint library.
 *
 *  Libint library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef INCLUDE_LIBINT2_LCAO_MOLDEN_H_
#define INCLUDE_LIBINT2_LCAO_MOLDEN_H_

#include <libint2/atom.h>
#include <libint2/basis.h>
#include <libint2/cgshell_ordering.h>
#include <libint2/chemistry/elements.h>
#include <libint2/shell.h>
#include <libint2/shgshell_ordering.h>

#include <cmath>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

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
  /// @tparam ShellSequence BasisSet or std::vector<Shell>
  /// @tparam Coeffs the type of LCAO coefficient matrix
  /// @tparam Energies the type of LCAO energy vector
  /// @tparam Occs the type of LCAO occupancy vector
  /// @param atoms the set of atoms (coordinates in atomic units)
  /// @param basis a sequence of shells; must meet Molden requirements (see
  /// below)
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
  /// @param bohr_to_angstrom the conversion factor from bohr to angstrom; the
  /// default is CODATA 2018 value
  /// @param coefficient_epsilon omit LCAO coefficients with absolute magnitude
  /// smaller than this value; set to 0 to write
  ///        all coefficients (some Molden parsers, e.g. Avogadro2, require
  ///        this)
  /// @throw std::logic_error if the basis does not conforms Molden
  ///        requirements
  /// @note Molden can only handle basis sets that:
  /// - p (l=1) shells are Cartesian, not solid harmonics
  /// - d, f, and g (l=2..4) shells are all Cartesian or all solid harmonics
  /// - there are no shells with l>5
  template <typename ShellSequence, typename Coeffs, typename Occs,
            typename Energies = Eigen::VectorXd>
  Export(
      const std::vector<Atom>& atoms, const ShellSequence& basis,
      const Coeffs& coefficients, const Occs& occupancies,
      const Energies& energies = Energies(),
      const std::vector<std::string>& symmetry_labels =
          std::vector<std::string>(),
      const std::vector<bool>& spincases = std::vector<bool>(),
      const double bohr_to_angstrom = constants::codata_2018::bohr_to_angstrom,
      double coefficient_epsilon = 5e-11)
      : atoms_(atoms),
        basis_(validate(basis)),
        coefs_(coefficients),
        occupancies_(occupancies),
        energies_(energies),
        labels_(symmetry_labels),
        spins_(spincases),
        bohr_to_angstrom_(bohr_to_angstrom),
        coefficient_epsilon_(coefficient_epsilon) {
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
      os << std::setw(4)
         << libint2::chemistry::get_element_info().at(Z - 1).symbol
         << std::setw(6) << (iatom + 1) << std::setw(6) << Z << std::setw(14)
         << atom.x << std::setw(14) << atom.y << std::setw(14) << atom.z
         << std::endl;
      ++iatom;
    }
  }

  /// writes the "[GTO]" section, as well as optional Cartesian/solid harmonics
  /// keywords, to ostream \c os
  void write_basis(std::ostream& os) const {
    bool f_found = false;
    os << "[GTO]" << std::endl;
    for (size_t iatom = 0; iatom < atoms_.size(); ++iatom) {
      os << std::setw(4) << (iatom + 1) << std::setw(4) << 0 << std::endl;
      for (auto ish : atom2shell_[iatom]) {
        const Shell& sh = basis_.at(ish);
        if (sh.contr.size() == 1) {
          const auto& contr = sh.contr[0];
          const auto l = contr.l;
          assert(l <= 4);  // only up to g functions are supported
          if (l == 3) f_found = true;
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
      // Molden default is cartesians throughout
      // dfg_is_cart_ is set to true even if there are no shells of a given type
      if (dfg_is_cart_[0]) {   // cartesian d
        if (!dfg_is_cart_[1])  // solid harmonic f
          os << "[7F]" << std::endl;
      } else {                   // solid harmonic d
        if (!dfg_is_cart_[1]) {  // solid harmonic f
          os << "[5D7F]" << std::endl;
        } else if (f_found) {  // cartesian f
          os << "[5D10F]" << std::endl;
        } else {  // no f functions
          os << "[5D]" << std::endl;
        }
      }
      if (!dfg_is_cart_[2])  // solid harmonic g
        os << "[9G]" << std::endl;
    }
  }

  /// writes the "[MO]" section to ostream \c os
  void write_lcao(std::ostream& os) const {
    os << "[MO]" << std::endl;
    for (int imo = 0; imo < coefs_.cols(); ++imo) {
      os << std::fixed << std::setprecision(10);
      os << std::setw(8) << "Sym= " << (labels_.empty() ? "A" : labels_.at(imo))
         << std::endl
         << std::setw(8) << "Ene= " << std::setw(16)
         << (energies_.rows() == 0 ? 0.0 : energies_(imo)) << std::endl
         << std::setw(8) << "Spin= "
         << (spins_.empty() ? "Alpha" : (spins_.at(imo) ? "Alpha" : "Beta"))
         << std::endl
         << std::setw(8) << "Occup= " << occupancies_(imo) << std::endl;
      os << std::scientific << std::uppercase << std::setprecision(10);
      for (int iao = 0; iao < coefs_.rows(); ++iao) {
        const auto C_ao_mo = coefs_(ao_map_[iao], imo);
        if (std::abs(C_ao_mo) >= coefficient_epsilon_) {
          os << std::setw(6) << (iao + 1) << " " << std::setw(16) << C_ao_mo
             << std::endl;
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

  double bohr_to_angstrom() const { return bohr_to_angstrom_; }

 private:
  const std::vector<Atom>& atoms_;
  const std::vector<Shell>& basis_;
  Eigen::MatrixXd coefs_;
  Eigen::VectorXd occupancies_;
  Eigen::VectorXd energies_;
  std::vector<std::string> labels_;
  std::vector<bool> spins_;
  double bohr_to_angstrom_;
  double coefficient_epsilon_;
  mutable bool dfg_is_cart_[3];  // whether {d, f, g} shells are cartesian
                                 // (true) or solid harmonics (false)
  std::vector<std::vector<long>>
      atom2shell_;  // maps atom -> shell indices in basis_
  std::vector<long>
      ao_map_;  // maps from the AOs ordered according to Molden
                // (atoms->shells, bf in shells ordered in the Molden order)
                // to the AOs ordered according to basis_

  /// @throw std::logic_error if the basis does not conforms Molden
  ///        requirements
  const std::vector<Shell>& validate(const std::vector<Shell>& shells) const {
    bool dfg_found[] = {false, false, false};
    for (int i = 0; i != sizeof(dfg_is_cart_) / sizeof(bool); ++i)
      dfg_is_cart_[i] = true;
    for (const auto& shell : shells) {
      for (const auto& contr : shell.contr) {
        if (contr.l > 4)
          throw std::logic_error(
              "molden::Export cannot handle shells with l > 4");

        switch (contr.l) {
          case 2:
          case 3:
          case 4: {
            if (!dfg_found[contr.l - 2]) {
              dfg_is_cart_[contr.l - 2] = !contr.pure;
              dfg_found[contr.l - 2] = true;
            }
            if (!contr.pure ^ dfg_is_cart_[contr.l - 2])
              throw std::logic_error(
                  "molden::Export only supports all-Cartesian or "
                  "all-solid-harmonics d/f/g shells");
          }

          default: {
          }  // l = 0 && 1 are fine
        }
      }
    }
    return shells;
  }

  /// @throw std::logic_error if the basis does not conforms Molden
  ///        requirements
  const std::vector<Shell>& validate(const BasisSet& bs) const {
    return validate(bs.shells());
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
            // Molden only handles Cartesian p shells, thus remap solid-harmonic
            // p shells to Cartesians:
            // - CCA: solid harmonics {y, z, x} -> Cartesian {x, y, z}
            // - Gaussian: solid harmonics {z, x, y} -> Cartesian {x, y, z}
            if (l == 1) {
#if LIBINT_SHGSHELL_ORDERING == LIBINT_SHGSHELL_ORDERING_STANDARD
              ao_map_[ao_molden] = ao + 2;
              ao_map_[ao_molden + 1] = ao;
              ao_map_[ao_molden + 2] = ao + 1;
#elif LIBINT_SHGSHELL_ORDERING == LIBINT_SHGSHELL_ORDERING_GAUSSIAN
              ao_map_[ao_molden] = ao + 1;
              ao_map_[ao_molden + 1] = ao + 2;
              ao_map_[ao_molden + 2] = ao;
#else
#error "unknown value of LIBINT_SHGSHELL_ORDERING"
#endif
              ao_molden += 3;
            } else {
              FOR_SOLIDHARM_MOLDEN(l, m)
              const auto ao_in_shell = libint2::INT_SOLIDHARMINDEX(l, m);
              ao_map_[ao_molden++] = ao + ao_in_shell;
              END_FOR_SOLIDHARM_MOLDEN
            }
            ao += 2 * l + 1;
          } else {
            long i, j, k;
            FOR_CART_MOLDEN(i, j, k, l)
            const auto ao_in_shell = libint2::INT_CARTINDEX(l, i, j);
            ao_map_[ao_molden++] = ao + ao_in_shell;
            END_FOR_CART_MOLDEN
            ao += INT_NCART(l);
          }
        }  // contraction loop
      }
    }
  }

};  // Export

/// Extension of the Molden exporter to support JMOL extensions for crystal
/// orbitals (see
/// <a>https://sourceforge.net/p/jmol/code/HEAD/tree/trunk/Jmol/src/org/jmol/adapter/readers/quantum/MoldenReader.java#l25</a>)
class PBCExport : public Export {
 public:
  /// @tparam Coeffs the type of LCAO coefficient matrix
  /// @tparam Energies the type of LCAO energy vector
  /// @tparam Occs the type of LCAO occupancy vector
  /// @param atoms the set of atoms (coordinates in atomic units)
  /// @param cell_axes the primitive vectors of the unit cell (in atomic units)
  /// @param basis the set of shells; must meet Molden requirements (see below)
  /// @param coefficients the matrix of LCAO coefficients (columns are LCAOs,
  ///        rows are AOs; AOs are ordered according to the order of shells in
  ///        \c basis and by the ordering conventions of this Libint
  ///        configuration)
  /// @param occupancies the vector of occupancies (size = # LCAOs)
  /// @param space_group (base-0) index of the space group in the International
  /// Tables of Crystallography (https://it.iucr.org/Ac/)
  /// @param energies the vector of energies (size = # of LCAOs); the default is
  ///        to assign zero to each LCAO
  /// @param symmetry_labels the vector of symmetry labels (size = # LCAOs); the
  ///        default is to assign empty label to each LCAO
  /// @param spincases the vector of spin cases (size = # LCAOs; true = spin-up
  ///        or m_s=1/2, false = spin-down or m_s=-1/2); the default is
  ///        to assign spin-up to each LCAO
  /// @param bohr_to_angstrom the conversion factor from bohr to angstrom; the
  /// default is CODATA 2018 value
  /// @throw std::logic_error if the basis does not conforms Molden
  ///        requirements
  /// @note Molden can only handle basis sets that:
  /// - p (l=1) shells are Cartesian, not solid harmonics
  /// - d, f, and g (l=2..4) shells are all Cartesian or all solid harmonics
  /// - there are no shells with l>5
  template <typename Coeffs, typename Occs, typename Energies = Eigen::VectorXd>
  PBCExport(
      const std::vector<Atom>& atoms,
      const std::array<Eigen::Vector3d, 3>& cell_axes,
      const std::vector<Shell>& basis, const Coeffs& coefficients,
      const Occs& occupancies, int space_group,
      const Energies& energies = Energies(),
      const std::vector<std::string>& symmetry_labels =
          std::vector<std::string>(),
      const std::vector<bool>& spincases = std::vector<bool>(),
      const double bohr_to_angstrom = constants::codata_2018::bohr_to_angstrom)
      : Export(atoms, basis, coefficients, occupancies, energies,
               symmetry_labels, spincases, bohr_to_angstrom),
        cell_axes_(cell_axes),
        space_group_(space_group) {
    // initialize_bf_map();
  }

  /// writes the "[SpaceGroup]" section to ostream \c os
  void write_space_group(std::ostream& os) const {
    os << "[SpaceGroup] (Number)" << std::endl;
    os << space_group_ << std::endl;
  }

  /// writes the "[Operators]" section to ostream \c os
  void write_operators(std::ostream& os) const {
    os << "[Operators]" << std::endl;
    os << "x, y, z" << std::endl;
  }

  /// writes the "[Cell]" section to ostream \c os
  void write_cell_axes(std::ostream& os) const {
    // https://sourceforge.net/p/jmol/code/HEAD/tree/trunk/Jmol/src/org/jmol/adapter/readers/quantum/MoldenReader.java#l107
    // suggests that [Cell] defaults to angstroms
    os << "[Cell]" << std::endl;
    {
      // convert vectors to abcɑβɣ
      const double a = cell_axes_[0].norm();
      const double b = cell_axes_[1].norm();
      const double c = cell_axes_[2].norm();
      const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
      const bool nonzero_a = a >= eps;
      const bool nonzero_b = b >= eps;
      const bool nonzero_c = c >= eps;
      static constexpr double right_angle = M_PI / 2;
      const double alpha =
          nonzero_b && nonzero_c
              ? std::acos(cell_axes_[1].dot(cell_axes_[2]) / (b * c))
              : right_angle;
      const double beta =
          nonzero_a && nonzero_c
              ? std::acos(cell_axes_[0].dot(cell_axes_[2]) / (a * c))
              : right_angle;
      const double gamma =
          nonzero_a && nonzero_b
              ? std::acos(cell_axes_[0].dot(cell_axes_[1]) / (a * b))
              : right_angle;
      const double radian_to_degree = 180 / M_PI;
      os << std::setw(12) << a * bohr_to_angstrom() << std::setw(12)
         << b * bohr_to_angstrom() << std::setw(12) << c * bohr_to_angstrom()
         << std::setw(12) << alpha * radian_to_degree << std::setw(12)
         << beta * radian_to_degree << std::setw(12) << gamma * radian_to_degree
         << std::endl;
    }
  }

  void write(const std::string& filename) const {
    std::ofstream os(filename);
    write_prologue(os);
    write_space_group(os);
    write_operators(os);
    write_cell_axes(os);
    write_atoms(os);
    write_basis(os);
    write_lcao(os);
  }

 private:
  std::array<Eigen::Vector3d, 3> cell_axes_;
  int space_group_;

};  // PBCExport

}  // namespace molden
}  // namespace libint2

#endif  // INCLUDE_LIBINT2_LCAO_MOLDEN_H_
