/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

#ifndef _libint2_src_bin_libint_utiltypes_h_
#define _libint2_src_bin_libint_utiltypes_h_

#include <string>

typedef enum {
  InBra=0, InKet=1
} FunctionPosition;
typedef enum {
  BraToKet=0, KetToBra=1
} FunctionMovement;
typedef enum {
  CartesianAxis_X=0, CartesianAxis_Y=1, CartesianAxis_Z=2
} CartesianAxis;

inline std::string to_string(CartesianAxis axis) {
  const char xyz_str[][2] = {"x", "y", "z"};
  return xyz_str[axis];
}

#endif /* header guard */
