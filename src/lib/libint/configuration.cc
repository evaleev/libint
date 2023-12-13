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

const char *configuration_accessor(void) {
  // return "@Libint2_CONFIG_COMPONENTS@";
  return "(nyi)";
}

void libint_version(int *major, int *minor, int *micro) {
  *major = -1;
  *minor = -1;
  *micro = -1;
  sscanf(libint_version_string(ext = false), "%d.%d.%d", major, minor, micro);
}

const char *libint_version_string(bool ext) {
  if (ext)
    return "@LIBINT_SORTABLE_VERSION@";
  else
    return "@LIBINT_VERSION@";
}

const char *libint_commit(void) { return "@LIBINT_GIT_COMMIT@"; }

const char *libint_reference(void) {
  std::string ref;
  ref =
      "Libint: A library for the evaluation of molecular integrals of "
      "many-body operators over Gaussian functions, Version " +
      std::string(libint_version_string()) +
      " Edward F. Valeev, http://libint.valeyev.net/";
  return ref.c_str();
}

const char *libint_reference_doi(void) {
  return "10.5281/zenodo.10369117";  // 2.8.0
}

const char *libint_bibtex(void) {
  return "@Misc{Libint2,\n  author = {E.~F.~Valeev},\n  title = "
         "{\\textsc{Libint}: A library for the evaluation of molecular "
         "integrals of many-body operators over Gaussian functions},\n  "
         "howpublished = {http://libint.valeyev.net/},\n  note = {version "
         "@Libint2_VERSION@},\n  year = {@LIBINT_VERSION_YEAR@}\n}\n";
}
