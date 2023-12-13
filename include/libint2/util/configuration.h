/*
 *  Copyright (C) 2004-2023 Edward F. Valeev
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

#ifndef _libint2_include_libint2_util_configuration_h_
#define _libint2_include_libint2_util_configuration_h_

/* Runtime accessor for the library configuration:
   integral derivatives, AM, orderings, etc.
   @return the semicolon-separated strings from CMake components */
const char* configuration_accessor();

#ifdef __cplusplus
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

namespace libint2 {

/// Runtime accessor for the library configuration:
/// integral derivatives, AM, orderings, etc.
/// @return the semicolon-separated strings from CMake components
inline std::string configuration_accessor() {
  std::string components = ::configuration_accessor();
  return components;
}

/// Runtime accessor for individual library configuration components:
/// integral derivatives, AM, orderings, etc.
/// @param[in] target CMake component with maximally uniform AM
/// @return whether target component available
inline bool supports(std::string component) {
  std::string segment;
  std::vector<std::string> seglist;
  std::stringstream ca(configuration_accessor());
  while (std::getline(ca, segment, ';')) {
    seglist.push_back(segment);
  }
  bool tf =
      (std::find(seglist.begin(), seglist.end(), component) != seglist.end());
  return tf;
}
}  // namespace libint2
#endif /* C++ guard */

#endif /* header guard */
