/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License, version 2,
 *  as published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
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

  /// reads the list of atoms from a file in the standard XYZ format supported by most chemistry software
  inline std::vector<Atom> read_dotxyz(std::istream& is) {
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
        std::cerr << "read_dotxyz: element symbol \"" << element_symbol << "\" is not recognized" << std::endl;
        throw "Did not recognize element symbol in .xyz file";
      }

      atoms[i].atomic_number = Z;

      // .xyz files report Cartesian coordinates in angstroms; convert to bohr
      const auto angstrom_to_bohr = 1 / 0.52917721092; // 2010 CODATA value
      //const auto angstrom_to_bohr = 1 / 0.529177249; // 1986 CODATA value, used by MPQC
      atoms[i].x = x * angstrom_to_bohr;
      atoms[i].y = y * angstrom_to_bohr;
      atoms[i].z = z * angstrom_to_bohr;
    }

    return atoms;
  }

  /// converts a vector of <code>Atom</code>s to a vector of point charges
  std::vector<std::pair<double, std::array<double, 3>>> make_point_charges(
      const std::vector<libint2::Atom>& atoms) {
    std::vector<std::pair<double, std::array<double, 3>>> q(atoms.size());
    for (const auto& atom : atoms) {
      q.emplace_back(static_cast<double>(atom.atomic_number),
                     std::array<double, 3>{{atom.x, atom.y, atom.z}});
    }
    return q;
  }

} // namespace libint2

#endif /* _libint2_src_lib_libint_atom_h_ */
