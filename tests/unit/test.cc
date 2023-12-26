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

#define CATCH_CONFIG_RUNNER

#include <libint2.hpp>

#include "catch.hpp"
#if !LIBINT2_CONSTEXPR_STATICS
#include <libint2/statics_definition.h>
#endif

int main(int argc, char* argv[]) {
  Catch::Session session;

  // add custom command line option to specify solid harmonics order

  std::string sho_str = "standard";

  // Build a new parser on top of Catch2's
  using namespace Catch::clara;
  auto cli =
      session.cli()  // Get Catch2's command line parser
      | Opt(sho_str,
            "shgshell-order")     // bind variable to a new option, with a hint
                                  // string
            ["--shgshell-order"]  // the option names it
                                  // will respond to
      ("solid harmonic order, valid values are \"standard\" (default) and "
       "\"gaussian\"");  // description string for the help output
  // Now pass the new composite back to Catch2 so it uses that
  session.cli(cli);

  // Let Catch2 (using Clara) parse the command line
  int returnCode = session.applyCommandLine(argc, argv);
  if (returnCode != 0)  // Indicates a command line error
    return returnCode;

  // global setup...
  libint2::SHGShellOrdering sho = libint2::SHGShellOrdering_Standard;
  if (sho_str == "gaussian") sho = libint2::SHGShellOrdering_Gaussian;
  libint2::set_solid_harmonics_ordering(sho);

  // initializes the Libint integrals library ... now ready to compute
  libint2::initialize();

  printf("Configuration S: sho=%d components=%s\n",
         libint2::solid_harmonics_ordering(),
         libint2::configuration_accessor().c_str());
  printf("Supports: dddd=%d mmmm=%d FF=%d\n", libint2::supports("eri_dddd_d0"),
         libint2::supports("eri_mmmm_d0"), libint2::supports("eri_FF_d0"));
  auto Mmp = libint2::libint_version();
  printf("Version: Numeric=%s Sortable=%s Commit=%s\n",
         libint2::libint_version_string(false).c_str(),
         libint2::libint_version_string(true).c_str(),
         libint2::libint_commit().c_str());
  printf("Version: Major=%d minor=%d patch=%d\n", std::get<0>(Mmp),
         std::get<1>(Mmp), std::get<2>(Mmp));
  printf("Citation: DOI=%s Ref=%s\n", libint2::libint_reference_doi().c_str(),
         libint2::libint_reference().c_str());
  printf("Citation: BibTeX=%s\n", libint2::libint_bibtex().c_str());

#ifdef LIBINT_HAS_MPFR
  // default to 256 bits of precision for mpf_class
  mpf_set_default_prec(256);
#endif

  int result = session.run();

  libint2::finalize();  // done with libint

  return result;
}
