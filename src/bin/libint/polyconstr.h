/*
 *  Copyright (C) 2004-2021 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_bin_libint_polyconstr_h_
#define _libint2_src_bin_libint_polyconstr_h_


namespace libint2 {

  /** ConstructablePolymorphically is a base for all objects
      which can be constructed using a SafePtr to a base or a
      SafePtr to ConstructablePolymorphically.
  */
  class ConstructablePolymorphically {
  protected:
    ConstructablePolymorphically() {}
    virtual ~ConstructablePolymorphically() {}
  };

};

#endif
