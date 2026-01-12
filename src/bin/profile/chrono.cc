/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint compiler.
 *
 *  Libint compiler is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint compiler is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint compiler.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <libint2/util/cxxstd.h>

#include <chrono>
#include <iomanip>
#include <iostream>

#if LIBINT2_CPLUSPLUS_STD < 2011
#error "compiler does not support C++11, try -std=c++11 flag"
#endif

int main(int argc, char* argv[]) {
  typedef std::chrono::high_resolution_clock clock_t;
  typedef std::chrono::time_point<clock_t> time_point_t;

  std::cout << "WARNING: turn off turboboost for reliable profiling"
            << std::endl;

  const size_t nrepeats = 20000000;
  std::cout << "nrepeats = " << nrepeats << std::endl;

  auto tstart = clock_t::now();
  for (size_t i = 0; i != nrepeats; ++i) {
    asm("nop");
  }
  auto tstop = clock_t::now();
  typedef std::chrono::duration<double> dur_t;
  const dur_t d0 = tstop - tstart;
  std::cout << "nop took " << std::setprecision(15)
            << d0.count() * 1e9 / nrepeats << " nanoseconds" << std::endl;

  tstart = clock_t::now();
  for (size_t i = 0; i != nrepeats; ++i) {
    const auto t = clock_t::now();
  }
  tstop = clock_t::now();
  const dur_t d1 = tstop - tstart;
  std::cout << "high_resolution_clock::now() took "
            << d1.count() * 1e9 / nrepeats << " nanoseconds" << std::endl;

  return 0;
}
