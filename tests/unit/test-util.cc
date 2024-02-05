/*
 *  Copyright (C) 2004-2024 Edward F. Valeev
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

#include <libint2/atom.h>

#include "catch.hpp"

using libint2::Atom;
using libint2::constants::codata_2014;
using libint2::constants::codata_2018;

TEST_CASE("XYZ reader", "[util]") {
  {  // fewer atoms than #atoms is not OK
    std::stringstream sstr;
    sstr << "2\n\nO 0 0 0\n";
    REQUIRE_THROWS_AS(libint2::read_dotxyz(sstr), std::logic_error);
  }
  {  // fewer atoms than #atoms+#axes is not OK
    std::stringstream sstr;
    sstr << "2\n\nO 0 0 0\nO 0 0 1\nAA 0 0 0\nBB 0 0 0\n";
    REQUIRE_THROWS_AS(libint2::read_dotxyz_pbc(sstr), std::logic_error);
  }
  {  // bad element symbol is not OK
    std::stringstream sstr;
    sstr << "2\n\nO 0 0 0\nZ 0 0 0\n";
    REQUIRE_THROWS_AS(libint2::read_dotxyz(sstr), std::logic_error);
  }
  {  // bad element symbol is not OK
    std::stringstream sstr;
    sstr << "2\n\n0 0 0 0\nO 0 0 0\n";
    REQUIRE_THROWS_AS(libint2::read_dotxyz(sstr), std::logic_error);
  }
  {  // duplicate cell parameters are not OK
    std::stringstream sstr;
    sstr << "2\n\nO 0 0 0\nO 0 0 1\nAA 2 0 0\nAA 0 2 0\nCC 0 0 2\n";
    REQUIRE_THROWS_AS(libint2::read_dotxyz_pbc(sstr), std::logic_error);
  }
  {  // missing cell parameters are not OK
    // not enough entries
    std::stringstream sstr1;
    sstr1 << "2\n\nO 0 0 0\nO 0 0 1\nAA 2 0 0\nCC 0 0 2\n";
    REQUIRE_THROWS_AS(libint2::read_dotxyz_pbc(sstr1), std::logic_error);
    // or too many atoms
    std::stringstream sstr2;
    sstr2 << "2\n\nO 0 0 0\nO 0 0 1\nO 0 0 3\nAA 2 0 0\nCC 0 0 2\n";
    REQUIRE_THROWS_AS(libint2::read_dotxyz_pbc(sstr2), std::logic_error);
  }
  {  // more atoms than #atoms is OK
    std::stringstream sstr;
    sstr << "2\n\nO 0 0 0\nO 0 0 0\nO 0 0 0\n";
    REQUIRE_NOTHROW(libint2::read_dotxyz(sstr));
  }
  {  // PBC input is OK as molecular input
    std::stringstream sstr;
    sstr << "2\n\nO 0 0 0\nO 0 0 1\nAA 2 0 0\nBB 0 2 0\nCC 0 0 2\n";
    REQUIRE_NOTHROW(libint2::read_dotxyz(sstr));
  }
  {  // validate results
    std::stringstream sstr;
    sstr << "2\n\nO 0 0 0\nO 0 0 1\n";
    auto atoms = libint2::read_dotxyz(sstr);
    const auto angstrom_to_bohr = 1 / codata_2018::bohr_to_angstrom;
    REQUIRE(atoms.size() == 2);
    REQUIRE(atoms[0] == Atom{8, 0., 0., 0.});
    REQUIRE(atoms[1] == Atom{8, 0., 0., 1. * angstrom_to_bohr});
  }
  {  // validate use of conversion factor
    std::stringstream sstr;
    sstr << "2\n\nO 0 0 0\nO 0 0 1\n";
    auto atoms = libint2::read_dotxyz(sstr, codata_2014::bohr_to_angstrom);
    const auto angstrom_to_bohr = 1 / codata_2014::bohr_to_angstrom;
    REQUIRE(atoms.size() == 2);
    REQUIRE(atoms[0] == Atom{8, 0., 0., 0.});
    REQUIRE(atoms[1] == Atom{8, 0., 0., 1. * angstrom_to_bohr});
  }
  {  // validate PBS results
    std::stringstream sstr;
    sstr
        << "2\n\nO 0 0 0\nO 0 0 1\nCC 0 0 4\nBB 0 3 0\nAA 2 0 0\n";  // axes can
                                                                     // be out
                                                                     // of order
    std::vector<Atom> atoms;
    std::array<std::array<double, 3>, 3> cell;
    std::tie(atoms, cell) = libint2::read_dotxyz_pbc(sstr);
    const auto angstrom_to_bohr = 1 / codata_2018::bohr_to_angstrom;
    REQUIRE(atoms.size() == 2);
    REQUIRE(atoms[0] == Atom{8, 0., 0., 0.});
    REQUIRE(atoms[1] == Atom{8, 0., 0., 1. * angstrom_to_bohr});
    REQUIRE(cell[0] == std::array<double, 3>{2. * angstrom_to_bohr, 0., 0.});
    REQUIRE(cell[1] == std::array<double, 3>{0., 3. * angstrom_to_bohr, 0.});
    REQUIRE(cell[2] == std::array<double, 3>{0., 0., 4. * angstrom_to_bohr});
  }
}
