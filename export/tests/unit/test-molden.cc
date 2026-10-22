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

#include <libint2/lcao/molden.h>

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "catch.hpp"
#include "fixture.h"

namespace {

/// @return the map from Molden AO index to the index of the AO in \p coefs,
///         as written by molden::Export::write_lcao for a single LCAO whose
///         coefficients are the (1-based) AO indices
std::vector<long> molden_ao_map(const std::vector<Atom>& atoms,
                                const std::vector<Shell>& shells) {
  const auto nao = libint2::nbf(shells);
  Eigen::MatrixXd coefs(nao, 1);
  for (long ao = 0; ao != nao; ++ao) coefs(ao, 0) = ao + 1;
  Eigen::VectorXd occs = Eigen::VectorXd::Ones(1);

  libint2::molden::Export xport(atoms, shells, coefs, occs);
  std::ostringstream oss;
  xport.write_lcao(oss);

  std::istringstream iss(oss.str());
  std::string line;
  // skip [MO], Sym=, Ene=, Spin=, Occup=
  for (int i = 0; i != 5; ++i) std::getline(iss, line);
  std::vector<long> result;
  long iao;
  double c;
  while (iss >> iao >> c) {
    REQUIRE(iao == static_cast<long>(result.size()) + 1);
    result.push_back(std::lround(c) - 1);
  }
  return result;
}

/// restores the solid harmonics ordering on scope exit
struct SHGShellOrderingGuard {
  libint2::SHGShellOrdering sho = libint2::solid_harmonics_ordering();
  ~SHGShellOrderingGuard() { libint2::set_solid_harmonics_ordering(sho); }
};

}  // namespace

TEST_CASE("Molden export follows runtime solid harmonics ordering",
          "[molden]") {
  const std::vector<Atom> atoms{{1, 0.0, 0.0, 0.0}};
  const std::vector<Shell> shells{
      Shell{{1.0}, {{1, true, {1.0}}}, {{0.0, 0.0, 0.0}}},  // solid harmonic p
      Shell{{1.0}, {{2, true, {1.0}}}, {{0.0, 0.0, 0.0}}}   // solid harmonic d
  };

  SHGShellOrderingGuard sho_guard;

  // Molden order: p = {x, y, z}, d = {m=0, +1, -1, +2, -2}
  SECTION("standard") {
    libint2::set_solid_harmonics_ordering(libint2::SHGShellOrdering_Standard);
    // libint: p = {y, z, x}, d = {m=-2, -1, 0, +1, +2}
    CHECK(molden_ao_map(atoms, shells) ==
          std::vector<long>{2, 0, 1, 5, 6, 4, 7, 3});
  }
  SECTION("gaussian") {
    libint2::set_solid_harmonics_ordering(libint2::SHGShellOrdering_Gaussian);
    // libint: p = {z, x, y}, d = {m=0, +1, -1, +2, -2}
    CHECK(molden_ao_map(atoms, shells) ==
          std::vector<long>{1, 2, 0, 3, 4, 5, 6, 7});
  }
}
