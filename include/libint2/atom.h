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

#ifndef _libint2_src_lib_libint_atom_h_
#define _libint2_src_lib_libint_atom_h_

#include <libint2/util/cxxstd.h>
#if LIBINT2_CPLUSPLUS_STD < 2011
#error "libint2/atom.h requires C++11 support"
#endif

#include <libint2/chemistry/elements.h>

#include <array>
#include <iostream>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>

namespace libint2 {

struct Atom {
  int atomic_number;
  double x, y, z;
};
inline bool operator==(const Atom &atom1, const Atom &atom2) {
  return atom1.atomic_number == atom2.atomic_number && atom1.x == atom2.x &&
         atom1.y == atom2.y && atom1.z == atom2.z;
}

namespace constants {
/// the 2022 CODATA reference set, available at
/// https://physics.nist.gov/cuu/pdf/wall_2022.pdf
struct codata_2022 {
  static constexpr double bohr_to_angstrom = 0.529177210544;
  static constexpr double angstrom_to_bohr = 1 / bohr_to_angstrom;
};
/// the 2018 CODATA reference set, available at
/// https://physics.nist.gov/cuu/pdf/wall_2018.pdf
struct codata_2018 {
  static constexpr double bohr_to_angstrom = 0.529177210903;
  static constexpr double angstrom_to_bohr = 1 / bohr_to_angstrom;
};
/// the 2014 CODATA reference set, available at DOI 10.1103/RevModPhys.88.035009
struct codata_2014 {
  static constexpr double bohr_to_angstrom = 0.52917721067;
  static constexpr double angstrom_to_bohr = 1 / bohr_to_angstrom;
};
/// the 2010 CODATA reference set, available at DOI 10.1103/RevModPhys.84.1527
struct codata_2010 {
  static constexpr double bohr_to_angstrom = 0.52917721092;
  static constexpr double angstrom_to_bohr = 1 / bohr_to_angstrom;
};
}  // namespace constants

}  // namespace libint2

namespace {

bool strcaseequal(const std::string &a, const std::string &b) {
  return a.size() == b.size() &&
         std::equal(a.begin(), a.end(), b.begin(), [](char a, char b) {
           return ::tolower(a) == ::tolower(b);
         });
}

/// reads the list of atoms from a file in the standard or PBC-extended XYZ
/// format \sa libint2::read_dotxyz \sa libint2::read_dotxyz_pbc
inline std::tuple<std::vector<libint2::Atom>,
                  std::array<std::array<double, 3>, 3>>
__libint2_read_dotxyz(std::istream &is, const double bohr_to_angstrom,
                      const bool pbc = false) {
  using libint2::Atom;
  const std::string caller =
      std::string("libint2::read_dotxyz") + (pbc ? "_pbc" : "");

  // first line = # of atoms
  size_t natom;
  is >> natom;

  // read off the rest of first line and discard
  std::string rest_of_line;
  std::getline(is, rest_of_line);

  // second line = comment
  std::string comment;
  std::getline(is, comment);

  // rest of lines are atoms (and unit cell parameters, if pbc = true)
  const auto nlines_expected = natom + (pbc ? 3 : 0);
  std::vector<Atom> atoms(natom, Atom{0, 0.0, 0.0, 0.0});
  std::array<std::array<double, 3>, 3> unit_cell({{{0.0, 0.0, 0.0}}});
  bool found_abc[3] = {false, false, false};
  for (size_t line = 0, atom_index = 0; line < nlines_expected; ++line) {
    if (is.eof())
      throw std::logic_error(caller + ": expected " +
                             std::to_string(nlines_expected) +
                             " sets of coordinates but only " +
                             std::to_string(line) + " received");

    // read line
    std::string linestr;
    std::getline(is, linestr);
    std::istringstream iss(linestr);
    // then parse ... this handles "extended" XYZ formats
    std::string element_symbol;
    double x, y, z;
    iss >> element_symbol >> x >> y >> z;

    // .xyz files report Cartesian coordinates in angstroms; convert to bohr
    const auto angstrom_to_bohr = 1 / bohr_to_angstrom;

    auto assign_atom = [angstrom_to_bohr](Atom &atom, int Z, double x, double y,
                                          double z) {
      atom.atomic_number = Z;
      atom.x = x * angstrom_to_bohr;
      atom.y = y * angstrom_to_bohr;
      atom.z = z * angstrom_to_bohr;
    };
    auto assign_xyz = [angstrom_to_bohr](std::array<double, 3> &xyz, double x,
                                         double y, double z) {
      xyz[0] = x * angstrom_to_bohr;
      xyz[1] = y * angstrom_to_bohr;
      xyz[2] = z * angstrom_to_bohr;
    };

    auto axis = -1;
    // if pbc = true, look for unit cell params
    if (pbc) {
      if (strcaseequal("AA", element_symbol)) axis = 0;
      if (strcaseequal("BB", element_symbol)) axis = 1;
      if (strcaseequal("CC", element_symbol)) axis = 2;
      if (axis != -1) {
        if (found_abc[axis])
          throw std::logic_error(
              caller + ": unit cell parameter along Cartesian axis " +
              std::to_string(axis) + " appears more than once");
        assign_xyz(unit_cell[axis], x, y, z);
        found_abc[axis] = true;
      }
    }

    // .xyz files report element labels, hence convert to atomic numbers
    if (axis == -1) {
      int Z = -1;
      for (const auto &e : libint2::chemistry::get_element_info()) {
        if (strcaseequal(e.symbol, element_symbol)) {
          Z = e.Z;
          break;
        }
      }
      if (Z == -1) {
        std::ostringstream oss;
        oss << caller << ": element symbol \"" << element_symbol
            << "\" is not recognized" << std::endl;
        throw std::logic_error(oss.str().c_str());
      }

      if (pbc &&
          atom_index == atoms.size()) {  // if PBC, check for too many atoms
        throw std::logic_error(caller + ": too many atoms");
      }
      assign_atom(atoms[atom_index++], Z, x, y, z);
    }
  }

  // make sure all 3 axes were specified
  if (pbc) {
    for (auto xyz = 0; xyz != 3; ++xyz)
      if (!found_abc[xyz]) {
        throw std::logic_error(caller +
                               ": unit cell parameter along Cartesian axis " +
                               std::to_string(xyz) + " not given");
      }
  }

  return std::make_tuple(atoms, unit_cell);
}

}  // anonymous namespace

namespace libint2 {

// clang-format off
/// reads the list of atoms from a file in the standard XYZ format supported
/// by most chemistry software (see <a href="https://en.wikipedia.org/wiki/XYZ_file_format">the Wikipedia page</a>)
/// \param[in] is the std::istream object from which the data will be read
/// \param[in] bohr_to_angstrom the conversion factor from Bohr (atomic unit
/// of length; Libint uses atomic units throughout) to angstrom (in which
/// the Cartesian coordinates are given in the XYZ file). The default is
/// the 2018 CODATA value given by the
/// libint2::constants::codata_2018::bohr_to_angstrom
/// constant.
/// \return a std::vector of Atom objects
/// \throw std::logic_error if cannot parse the contents of \c is
// clang-format off
inline std::vector<Atom> read_dotxyz(
    std::istream &is,
    const double bohr_to_angstrom = constants::codata_2018::bohr_to_angstrom) {
  std::vector<Atom> atoms;
  std::tie(atoms, std::ignore) =
      __libint2_read_dotxyz(is, bohr_to_angstrom, false);
  return atoms;
}

/// reads the list of atoms from a file in the PBC-extended XYZ format
/// \note The unit cell vectors in PBC-extended XYZ file are specified as atoms with
///       element symbols "AA", "BB", and "CC" (N.B. the element symbols are not
///       case-sensitive).
/// Omitting all three unit cell vectors is equivalent to an infinite unit cell (no periodicity in any direction).
/// \param[in] is the std::istream object from which the data will be read
/// \param[in] bohr_to_angstrom the conversion factor from Bohr (atomic unit
/// of length; Libint uses atomic units throughout) to angstrom (in which
/// the Cartesian coordinates are given in the XYZ file). The default is
/// the 2018 CODATA value given by the
/// libint2::constants::codata_2018::bohr_to_angstrom
/// constant.
/// \return a tuple composed of the list of atoms and an array of 3
///         unit cell vectors, \c A , \c B , and \c C .
/// \throw std::logic_error if cannot parse the contents of \c is
inline auto read_dotxyz_pbc(
    std::istream &is,
    const double bohr_to_angstrom = constants::codata_2018::bohr_to_angstrom)
    -> decltype(__libint2_read_dotxyz(is, bohr_to_angstrom, true)) {
  return __libint2_read_dotxyz(is, bohr_to_angstrom, true);
}

/// converts a vector of <code>Atom</code>s to a vector of point charges
std::vector<std::pair<double, std::array<double, 3>>> inline make_point_charges(
    const std::vector<libint2::Atom> &atoms) {
  std::vector<std::pair<double, std::array<double, 3>>> q;
  q.reserve(atoms.size());
  for (const auto &atom : atoms) {
    q.emplace_back(static_cast<double>(atom.atomic_number),
                   std::array<double, 3>{{atom.x, atom.y, atom.z}});
  }
  return q;
}

} // namespace libint2

#endif /* _libint2_src_lib_libint_atom_h_ */
