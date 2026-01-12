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

#include <iter.h>

using namespace std;
using namespace libint2;

//
// SubIterator
//
SubIterator::SubIterator() {}

SubIterator::~SubIterator() {}

const ConstructablePolymorphically& SubIterator::pelem() const {
  throw ProgrammingError(
      "virtual ConstructablePolymorphically& SubIterator::pelem() -- used but "
      "not overloaded");
}

#ifdef ALLOW_PARTIALLY_SPECIALIZED_NESTED_TEMPLATES

//
// SubIteratorBase::init_subobj<CGF>
//
template <class T, class P>
template <>
void SubIteratorBase<T, P>::init_subobj<CGF>() {
  subobj_.push_back(obj_);
}

template <class T, class P>
template <>
void SubIteratorBase<T, P>::delete_subobj<CGF>() {}

//
// SubIteratorBase::init_subobj<CGShell>
//
template <class T, class P>
template <>
void SubIteratorBase<T, P>::init_subobj<CGShell>() {
  P::cgshell_to_cgfvector(obj_, subobj_);
}

template <class T, class P>
template <>
void SubIteratorBase<T, P>::delete_subobj<CGShell>() {
  int nelem = subobj_.size();
  for (int i = 0; i < nelem; i++) subobj_[i]->~CGF();
}

#endif
