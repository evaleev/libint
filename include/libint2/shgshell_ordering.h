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

#ifndef _libint2_src_bin_libint_shgshellordering_h_
#define _libint2_src_bin_libint_shgshellordering_h_

#include <libint2/config.h>

#include <cmath>

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

/* Computes an index to a Cartesian function within a shell given
 * l = total angular momentum
 * m = real solid harmonic index (|m| = the absolute value of the projection of
 * the angular momentum on the z axis) m runs from -l to l
 */
namespace libint2 {
inline int INT_SOLIDHARMINDEX_STANDARD(int l, int m) { return m + l; }
}

/* This sets up the above loop over cartesian exponents as follows
 * int m;
 * FOR_SOLIDHARM_STANDARD(l,m)
 * END_FOR_SOLIDHARM
 */
#define FOR_SOLIDHARM_STANDARD(l, m) for ((m) = -(l); (m) <= (l); ++(m)) {
#define END_FOR_SOLIDHARM }

/* Computes an index to a Cartesian function within a shell given
 * l = total angular momentum
 * m = real solid harmonic index (|m| = the absolute value of the projection of
 * the angular momentum on the z axis) m runs as 0, +1, -1, +2, -2 ... +l, -l
 */
namespace libint2 {
inline int INT_SOLIDHARMINDEX_GAUSSIAN(int l, int m) {
  return 2 * std::abs(m) + (m > 0 ? -1 : 0);
}
}

/* This sets up the above loop over cartesian exponents as follows
 * int m;
 * FOR_SOLIDHARM_GAUSSIAN(l,m)
 * END_FOR_SOLIDHARM
 */
#define FOR_SOLIDHARM_GAUSSIAN(l, m) \
  for ((m) = 0; (m) != (l) + 1; (m) = ((m) > 0 ? -(m) : 1 - (m))) {
/// these always-available macros encode orderings assumed by Molden ... useful
/// for writing Molden-specific code
/// @note Molden and Gaussian ordering coincide

namespace libint2 {
inline int INT_SOLIDHARMINDEX_MOLDEN(int l, int m) {
  return INT_SOLIDHARMINDEX_GAUSSIAN(l, m);
}
}

#define FOR_SOLIDHARM_MOLDEN(l, m) FOR_SOLIDHARM_GAUSSIAN(l, m)
#define END_FOR_SOLIDHARM_MOLDEN END_FOR_SOLIDHARM

namespace libint2 {
LIBINT_DEPRECATED(
    "please use "
    "libint2::INT_SOLIDHARMINDEX(libint2::solid_harmonics_ordering(), l, m) "
    "instead. Current function returns the standard or gaussian index based on "
    "build configuration, not on libint2::set_solid_harmonics_ordering() "
    "argument.")
inline int INT_SOLIDHARMINDEX(int l, int m) {
#if LIBINT_SHGSHELL_ORDERING == LIBINT_SHGSHELL_ORDERING_STANDARD
  return libint2::INT_SOLIDHARMINDEX_STANDARD(l, m);
#elif LIBINT_SHGSHELL_ORDERING == LIBINT_SHGSHELL_ORDERING_GAUSSIAN
  return libint2::INT_SOLIDHARMINDEX_GAUSSIAN(l, m);
#else
#error "unknown value of macro LIBINT_SHGSHELL_ORDERING"
#endif
}

inline int INT_SOLIDHARMINDEX(int sho, int l, int m) {
  if (sho == LIBINT_SHGSHELL_ORDERING_STANDARD)
    return libint2::INT_SOLIDHARMINDEX_STANDARD(l, m);
  else
    return libint2::INT_SOLIDHARMINDEX_GAUSSIAN(l, m);
}
}

// FOR_SOLIDHARM dispatches to the appropriate FOR_SOLIDHARM_ macro
#if LIBINT_SHGSHELL_ORDERING == LIBINT_SHGSHELL_ORDERING_STANDARD
#define FOR_SOLIDHARM FOR_SOLIDHARM_STANDARD
#elif LIBINT_SHGSHELL_ORDERING == LIBINT_SHGSHELL_ORDERING_GAUSSIAN
#define FOR_SOLIDHARM FOR_SOLIDHARM_GAUSSIAN
#else
#error "unknown value of macro LIBINT_SHGSHELL_ORDERING"
#endif

#endif  // _libint2_src_bin_libint_shgshellordering_h_
