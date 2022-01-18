/*
 *  Copyright (C) 2018-2021 Edward F. Valeev
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

#include "../unit/catch.hpp"
#include <libint2.hpp>

#ifdef FC_DUMMY_MAIN
#  ifdef __cplusplus
extern "C"
#  endif
int FC_DUMMY_MAIN() { return 1; }
#endif

int main( int argc, char* argv[] )
{
  Catch::Session session;
  // global setup...
  // initializes the Libint integrals library ... now ready to compute
  libint2::initialize();

  int result = session.run( argc, argv );

  libint2::finalize(); // done with libint

  return result;
}
