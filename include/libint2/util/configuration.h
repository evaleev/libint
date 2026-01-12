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

#ifndef _libint2_include_libint2_util_configuration_h_
#define _libint2_include_libint2_util_configuration_h_

#include <stdbool.h>

/* Runtime accessor for the library configuration:
   integral derivatives, AM, orderings, etc.
   @return the semicolon-separated strings from CMake components */
const char *configuration_accessor();

/* Get the major, minor, and micro version of Libint */
void libint_version(int *, int *, int *);

/* Get the version of Libint as a string
    @return the version string. At release, strictly "M.m.p" (no alpha/rc/etc.).
    Beyond release (arg=true), returns "M.m.p.postD" where D is distance from
   release. Beyond release (arg=false), returns most recent release, "M.m.p". */
const char *libint_version_string(bool);

/* Get the git commit at which library was generated
    @return the commit as a 7-char abbreviated string */
const char *libint_commit(void);

/* Literature citation
    @return the citation string including description and version */
const char *libint_reference(void);

/* Literature citation DOI
   @return the string of DOI for latest tag */
const char *libint_reference_doi(void);

/* BibTeX for citing Libint
   @return the string for literature citation */
const char *libint_bibtex(void);

#ifdef __cplusplus
#include <algorithm>
#include <sstream>
#include <string>
#include <tuple>
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

/// Get the major, minor, and micro version of Libint */
/// @return the components of the last release
inline std::tuple<int, int, int> libint_version() {
  int vmajor, vminor, vmicro;
  ::libint_version(&vmajor, &vminor, &vmicro);
  return std::make_tuple(vmajor, vminor, vmicro);
}

/// Get the version of Libint as a string
/// @param[in] whether to return the simple-sortable last release or a
/// per-commit version
/// @return the version string. At release, strictly "M.m.p" (no alpha/rc/etc.).
/// Beyond release (ext=true), returns "M.m.p.postD" where D is distance from
/// release. Beyond release (ext=false), returns most recent release, "M.m.p".
inline std::string libint_version_string(bool ext = true) {
  std::string version = ::libint_version_string(ext);
  return version;
}

/// Get the git commit at which library was generated
/// @return the commit as a 7-char abbreviated string
inline std::string libint_commit(void) {
  std::string commit = ::libint_commit();
  return commit;
}

/// Literature citation
/// @return the citation string including description and version
inline std::string libint_reference(void) {
  std::string ref = ::libint_reference();
  return ref;
}

/// Literature citation DOI
/// @return the string of DOI for latest tag
inline std::string libint_reference_doi(void) {
  std::string ref = ::libint_reference_doi();
  return ref;
}

/// BibTeX for citing Libint
/// @return the string for literature citation
inline std::string libint_bibtex(void) {
  std::string ref = ::libint_bibtex();
  return ref;
}

}  // namespace libint2
#endif /* C++ guard */

#endif /* header guard */
