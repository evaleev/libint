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

#ifndef LIBINT_FIXTURE_H
#define LIBINT_FIXTURE_H

#include <libint2.hpp>

using libint2::Atom;
using libint2::BasisSet;
using libint2::BraKet;
using libint2::CartesianShellNormalization;
using libint2::Engine;
using libint2::Operator;
using libint2::Shell;

namespace libint2 {
namespace unit {

class DefaultFixture {
 public:
  DefaultFixture()
      : atoms{{8, 0., 0., 0.},
              {8, 0., 0., 2.},
              {1, 0., -1., -1.},
              {1, 0., 1., 3.}},
        obs("6-31g*", atoms),
        dfbs("aug-cc-pvdz", atoms) {}

 protected:
  std::vector<Atom> atoms;
  BasisSet obs, dfbs;
};

}  // namespace unit
}  // namespace libint2

#endif  // LIBINT_FIXTURE_H
