/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

#include <iostream>
#include <libint/util.h>
#include <libint/intpart_iter.h>
#include <libint/buildtest.h>
#include <libint/bfset.h>
#include <libint/master_ints_list.h>

using namespace libint2;

long living_count = 0;

namespace {
  int try_main (int argc, char* argv[]);
};

int main (int argc, char* argv[])
{
  try {
    try_main(argc,argv);
  }
  catch(std::exception& a) {
    std::cerr << std::endl
         << "  WARNING! Caught a standard exception:" << std::endl
         << "    " << a.what() << std::endl << std::endl;
  }
  return 0;
}

namespace {

  unsigned int am[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  typedef TwoPRep_11_11_sq ERIQtet;

  int try_main (int argc, char* argv[])
  {
#if LIBINT_CONTRACTED_INTS
    CGShell::set_contracted_default_value(true);
#else
    CGShell::set_contracted_default_value(false);
#endif
    TesterCmdLine<4> cmdline(argc,argv);
    const std::vector<unsigned int>& l = cmdline.am();

    // set this to the order of deriv ERI
    // will pick up from command line eventually
    const unsigned int deriv_order = 0;
    CartesianDerivIterator<4u> diter(deriv_order);

    //
    // generate all targets
    //
    std::vector< SafePtr<ERIQtet> > targets;
    // iterate over the distribution of all derivative quanta among 12 different directions
    bool last_deriv = false;
    do {
      CGShell a(l[0]);  for(unsigned int xyz=0; xyz<3; ++xyz) a.deriv().inc(xyz, diter.value(0+xyz));
      CGShell b(l[1]);  for(unsigned int xyz=0; xyz<3; ++xyz) b.deriv().inc(xyz, diter.value(3+xyz));
      CGShell c(l[2]);  for(unsigned int xyz=0; xyz<3; ++xyz) c.deriv().inc(xyz, diter.value(6+xyz));
      CGShell d(l[3]);  for(unsigned int xyz=0; xyz<3; ++xyz) d.deriv().inc(xyz, diter.value(9+xyz));

      const SafePtr<ERIQtet> abcd0 = ERIQtet::Instance(a,b,c,d, 0u);
      targets.push_back(abcd0);
      last_deriv = diter.last();
      if (!last_deriv) diter.next();
    } while (!last_deriv);

    BuildTest<ERIQtet,true>(targets,
			    cmdline.size_to_unroll(),
			    cmdline.veclen(),
			    cmdline.vectorize_by_line(),
			    cmdline.do_cse(),
			    "eri0",
			    std::cout);

    return 0;
  }
};

