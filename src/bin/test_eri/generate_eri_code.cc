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
#include <libint/master_ints_list.h>
#include <libint/util.h>
#include <libint2/deriv_iter.h>

#include <iostream>

using namespace libint2;

long living_count = 0;

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

typedef TwoPRep_11_11_sq ERIQtet;

int try_main(int argc, char* argv[]) {
#if LIBINT_CONTRACTED_INTS
  CGShell::set_contracted_default_value(true);
#else
  CGShell::set_contracted_default_value(false);
#endif
  TesterCmdLine<4> cmdline(argc, argv);
  const std::vector<unsigned int>& l = cmdline.am();

  // set this to the order of deriv ERI
  // will pick up from command line eventually
  const unsigned int deriv_order = 0;
  CartesianDerivIterator<4u> diter(deriv_order);

  //
  // generate all targets
  //
  std::vector<std::shared_ptr<ERIQtet> > targets;
  // iterate over the distribution of all derivative quanta among 12 different
  // directions
  bool last_deriv = false;
  do {
    CGShell a(l[0]);
    for (unsigned int xyz = 0; xyz < 3; ++xyz)
      a.deriv().inc(xyz, (*diter)[0 + xyz]);
    CGShell b(l[1]);
    for (unsigned int xyz = 0; xyz < 3; ++xyz)
      b.deriv().inc(xyz, (*diter)[3 + xyz]);
    CGShell c(l[2]);
    for (unsigned int xyz = 0; xyz < 3; ++xyz)
      c.deriv().inc(xyz, (*diter)[6 + xyz]);
    CGShell d(l[3]);
    for (unsigned int xyz = 0; xyz < 3; ++xyz)
      d.deriv().inc(xyz, (*diter)[9 + xyz]);

    const std::shared_ptr<ERIQtet> abcd0 = ERIQtet::Instance(a, b, c, d, 0u);
    targets.push_back(abcd0);
    last_deriv = diter.last();
    if (!last_deriv) diter.next();
  } while (!last_deriv);

  BuildTest<ERIQtet, true>(targets, cmdline.size_to_unroll(), cmdline.veclen(),
                           cmdline.vectorize_by_line(), cmdline.do_cse(),
                           "eri0", std::cout);

  return 0;
}
};  // namespace
