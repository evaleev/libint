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

#include <global_macros.h>
#include <libint2/cgshell_ordering.h>
#include <policy.h>
#include <rr.h>

using namespace std;

/*
 Definition of a generic StdLibintTDPolicy is provided in policy.h
 */

namespace libint2 {

/**
StdLibintTDPolicy<CGShell>::init_subobj initializes CGFs in canonical order.
 */

template <>
void StdLibintTDPolicy<CGShell>::init_subobj(
    const StdLibintTDPolicy<CGShell>::obj_stype& cgshell,
    vector<StdLibintTDPolicy<CGShell>::subobj_stype>& cgfs) {
  if (cgshell.is_unit()) {
    cgfs.push_back(CGF::unit());
  } else {
    unsigned int am = TypeTraits<CGShell>::const_ref(cgshell).qn();
    unsigned int qn[3] = {0, 0, 0};
    int lx, ly, lz;
    FOR_CART(lx, ly, lz, am)
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
void StdLibintTDPolicy<CGShell>::dealloc_subobj(
    vector<StdLibintTDPolicy<CGShell>::subobj_stype>& subobj) {}

};  // namespace libint2
