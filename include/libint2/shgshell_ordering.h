/*
 *  Copyright (C) 2004-2020 Edward F. Valeev
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

#ifndef _libint2_src_bin_libint_shgshellordering_h_
#define _libint2_src_bin_libint_shgshellordering_h_

#include <cmath>

#include <libint2/config.h>

namespace libint2 {

enum SHGShellOrdering {
  SHGShellOrdering_Standard = LIBINT_SHGSHELL_ORDERING_STANDARD,
  SHGShellOrdering_Gaussian = LIBINT_SHGSHELL_ORDERING_GAUSSIAN,
  SHGShellOrdering_MOLDEN  // same as Gaussian
};

}

//
// Macros that define orderings
//

#if LIBINT_SHGSHELL_ORDERING == LIBINT_SHGSHELL_ORDERING_STANDARD

/* Computes an index to a Cartesian function within a shell given
 * l = total angular momentum
 * m = real solid harmonic index (|m| = the absolute value of the projection of
 * the angular momentum on the z axis) m runs from -l to l
 */
#define INT_SOLIDHARMINDEX(l, m) ((m) + (l))

/* This sets up the above loop over cartesian exponents as follows
 * int m;
 * FOR_SOLIDHARM(l,m)
 * END_FOR_SOLIDHARM
 */
#define FOR_SOLIDHARM(l, m) for ((m) = -(l); (m) <= (l); ++(m)) {
#define END_FOR_SOLIDHARM }

#endif  // Standard ordering

#if LIBINT_SHGSHELL_ORDERING == LIBINT_SHGSHELL_ORDERING_GAUSSIAN

/* Computes an index to a Cartesian function within a shell given
 * l = total angular momentum
 * m = real solid harmonic index (|m| = the absolute value of the projection of
 * the angular momentum on the z axis) m runs as 0, +1, -1, +2, -2 ... +l, -l
 */
#define INT_SOLIDHARMINDEX(l, m) (2 * std::abs(m) + ((m) > 0 ? -1 : 0))

/* This sets up the above loop over cartesian exponents as follows
 * int m;
 * FOR_SOLIDHARM(l,m)
 * END_FOR_SOLIDHARM
 */
#define FOR_SOLIDHARM(l, m) \
  for ((m) = 0; (m) != (l) + 1; (m) = ((m) > 0 ? -(m) : 1 - (m))) {
#define END_FOR_SOLIDHARM }

#endif  // Gaussian ordering

/// these always-available macros encode orderings assumed by Molden

#define INT_SOLIDHARMINDEX_MOLDEN(l, m) (2 * std::abs(m) + ((m) > 0 ? -1 : 0))
#define FOR_SOLIDHARM_MOLDEN(l, m) \
  for ((m) = 0; (m) != (l) + 1; (m) = ((m) > 0 ? -(m) : 1 - (m))) {
#define END_FOR_SOLIDHARM_MOLDEN }

#endif  // _libint2_src_bin_libint_shgshellordering_h_
