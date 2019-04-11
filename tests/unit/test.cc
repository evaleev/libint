/*
 *  Copyright (C) 2018 Edward F. Valeev
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

#define CATCH_CONFIG_RUNNER

#include "catch.hpp"
#include <libint2.hpp>
#if !LIBINT2_CONSTEXPR_STATICS
#  include <libint2/statics_definition.h>
#endif

int main( int argc, char* argv[] )
{
  Catch::Session session;
  // global setup...
  // initializes the Libint integrals library ... now ready to compute
  libint2::initialize();

#ifdef LIBINT_HAS_MPFR
  // default to 256 bits of precision for mpf_class
  mpf_set_default_prec(256);
#endif

  int result = session.run( argc, argv );

  libint2::finalize(); // done with libint

  return result;
}
