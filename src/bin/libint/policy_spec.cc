/*
 *  Copyright (C) 2004-2019 Edward F. Valeev
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

#include <policy.h>
#include <rr.h>
#include <libint2/cgshell_ordering.h>

using namespace std;

/*
 Definition of a generic StdLibintTDPolicy is provided in policy.h
 */

namespace {
  /// returns an xyz which is not a and not b
  int notxyz(int a, int b);
  /// returns an ordered pair of xyz (i.e. xy, not yx) which is not a
  std::pair<int,int> notxyz(int a);
}

namespace libint2 {

  /**
  StdLibintTDPolicy<CGShell>::init_subobj initializes CGFs in canonical order.
   */

template <>
void
StdLibintTDPolicy<CGShell>::init_subobj(const StdLibintTDPolicy<CGShell>::obj_stype& cgshell,
                                        vector<StdLibintTDPolicy<CGShell>::subobj_stype>& cgfs)
{
  if (cgshell.is_unit()) {
      cgfs.push_back(CGF::unit());
  }
  else {
    unsigned int am = TypeTraits<CGShell>::const_ref(cgshell).qn();
    unsigned int qn[3] = {0, 0, 0};
    int lx, ly, lz;
    FOR_CART(lx,ly,lz,am)
      qn[0] = lx;
      qn[1] = ly;
      qn[2] = lz;
      subobj_stype cgf(qn, cgshell.pure_sh());
      cgf.deriv() = cgshell.deriv();
      if (cgshell.contracted()) cgf.contract();
      cgfs.push_back(cgf);
    END_FOR_CART
  }
}

template <>
void
StdLibintTDPolicy<CGShell>::dealloc_subobj(vector<StdLibintTDPolicy<CGShell>::subobj_stype>& subobj)
{
}

};

////

namespace {
  int notxyz(int a, int b) {
    if (a == b)
      throw libint2::ProgrammingError("notxyz(a,b) -- a equals b");
    int amax = std::max(a,b);
    int amin = std::min(a,b);
    if (amin == 0 && amax == 1)
      return 2;
    if (amin == 0 && amax == 2)
      return 1;
    if (amin == 1 && amax == 2)
      return 0;
    abort(); // unreachable
  }

  std::pair<int,int> notxyz(int a) {
    switch(a) {
    case 0: return make_pair(1,2); break;
    case 1: return make_pair(0,2); break;
    case 2: return make_pair(0,1); break;
    }
    abort(); // unreachable
  }
}
