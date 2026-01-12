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

#include <purgeable.h>

#include <algorithm>
#include <cassert>

using namespace libint2;

PurgeableStacks* PurgeableStacks::instance_ = 0;

PurgeableStacks* PurgeableStacks::Instance() {
  if (instance_ == 0) instance_ = new PurgeableStacks;
  return instance_;
}

void PurgeableStacks::register_stack(stack_type* stack) {
  typedef std::vector<stack_type*>::const_iterator citer;
  citer v = std::find(stacks_.begin(), stacks_.end(), stack);
  assert(v == stacks_.end());
  stacks_.push_back(stack);
}

void PurgeableStacks::purge() {
  typedef std::vector<stack_type*>::iterator iter;
  for (iter s = stacks_.begin(); s != stacks_.end(); ++s) (*s)->purge();
}
