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

#include <libint/bfset.h>
#include <libint/buildtest.h>

#include <iostream>

using namespace libint2;

namespace {
int try_main(int argc, char* argv[]);
};

int main(int argc, char* argv[]) {
  try {
    try_main(argc, argv);
  } catch (std::exception& a) {
    std::cerr << std::endl
              << "  WARNING! Caught a standard exception:" << std::endl
              << "    " << a.what() << std::endl
              << std::endl;
  }
  return 0;
}

namespace {

unsigned int am[][1] = {{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}};
typedef TiG12_11_11<CGShell, 0> IntegralSet;

int try_main(int argc, char* argv[]) {
  TesterCmdLine<4> cmdline(argc, argv);
  std::vector<unsigned int> l = cmdline.am();
  const std::shared_ptr<IntegralSet> integral =
      IntegralSet::Instance(CGShell(am[l[0]]), CGShell(am[l[1]]),
                            CGShell(am[l[2]]), CGShell(am[l[3]]));
  BuildTest<IntegralSet, true>(integral, cmdline.size_to_unroll(),
                               cmdline.veclen(), cmdline.vectorize_by_line(),
                               cmdline.do_cse(), "t0g12", std::cout);
  return 0;
}
};  // namespace
