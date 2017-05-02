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

#ifndef _libint2_src_lib_libint_atom_h_
#define _libint2_src_lib_libint_atom_h_

#include <libint2/util/cxxstd.h>
#if LIBINT2_CPLUSPLUS_STD < 2011
# error "libint2/atom.h requires C++11 support"
#endif

#include <array>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

#include <libint2/chemistry/elements.h>

namespace libint2 {

  struct Atom {
      int atomic_number;
      double x, y, z;
  };

  namespace {
    bool strcaseequal(const std::string& a, const std::string& b) {
      return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin(),
                                                [](char a, char b) {return ::tolower(a) == ::tolower(b);}
                                               );
    }
  }

  namespace constants {
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

  /// reads the list of atoms from a file in the standard XYZ format supported
  /// by most chemistry software (see <a
  /// href="https://en.wikipedia.org/wiki/XYZ_file_format">the Wikipedia
  /// page</a>)
  /// \param is[in] the std::istream object from which the data will be read
  /// \param bohr_to_angstrom[in] the conversion factor from Bohr (atomic unit
  /// of length; Libint uses atomic units throughout) to angstrom (in which
  /// the Cartesian coordinates are given in the XYZ file). The default is
  /// the 2010 CODATA value given by the
  /// libint2::constants::codata_2010::bohr_to_angstrom
  /// constant.
  /// \return a std::vector of Atom objects
  /// \throw std::runtime_error if cannot parse the contents of \c is
  inline std::vector<Atom> read_dotxyz(
      std::istream& is,
      const double bohr_to_angstrom = constants::codata_2010::bohr_to_angstrom) {
    const double angstrom_to_bohr = 1 / bohr_to_angstrom;
    // first line = # of atoms
    size_t natom;
    is >> natom;
    // read off the rest of first line and discard
    std::string rest_of_line;
    std::getline(is, rest_of_line);

    // second line = comment
    std::string comment;
    std::getline(is, comment);

    // rest of lines are atoms
    std::vector<Atom> atoms(natom);
    for (auto i = 0; i < natom; i++) {
      // read line
      std::string line;
      std::getline(is, line);
      std::istringstream iss(line);
      // then parse ... this handles "extended" XYZ formats
      std::string element_symbol;
      double x, y, z;
      iss >> element_symbol >> x >> y >> z;

      // .xyz files report element labels, hence convert to atomic numbers
      int Z = -1;
      using libint2::chemistry::element_info;
      for(const auto& e: element_info) {
        if (strcaseequal(e.symbol, element_symbol)) {
          Z = e.Z;
          break;
        }
      }
      if (Z == -1) {
        std::ostringstream oss;
        oss << "read_dotxyz: element symbol \"" << element_symbol << "\" is not recognized" << std::endl;
        throw std::runtime_error(oss.str().c_str());
      }

      atoms[i].atomic_number = Z;

      // .xyz files report Cartesian coordinates in angstroms; convert to bohr
      atoms[i].x = x * angstrom_to_bohr;
      atoms[i].y = y * angstrom_to_bohr;
      atoms[i].z = z * angstrom_to_bohr;
    }

    return atoms;
  }

  /// converts a vector of <code>Atom</code>s to a vector of point charges
  std::vector<std::pair<
      double,
      std::array<double, 3>>> inline make_point_charges(const std::
                                                            vector<
                                                                libint2::Atom>&
                                                                atoms) {
    std::vector<std::pair<double, std::array<double, 3>>> q(atoms.size());
    for (const auto& atom : atoms) {
      q.emplace_back(static_cast<double>(atom.atomic_number),
                     std::array<double, 3>{{atom.x, atom.y, atom.z}});
    }
    return q;
  }

} // namespace libint2

#endif /* _libint2_src_lib_libint_atom_h_ */
