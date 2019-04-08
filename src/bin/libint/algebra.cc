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

#include <algebra.h>
#include <global_macros.h>

namespace libint2 {

  template <>
  void
  AlgebraicOperator<DGVertex>::add_exit_arc(const SafePtr<DGArc>& a)
  {
    DGVertex::add_exit_arc(a);
#if CHECK_SAFETY
    if (num_exit_arcs() > 2)
      throw std::runtime_error("AlgebraicOperator<DGVertex>::add_exit_arc() -- number of exit arcs is now greater than 2!");
#endif
    if (left_->equiv(a->dest()))
      left_ = a->dest();
    else if (right_->equiv(a->dest()))
      right_ = a->dest();
    else
      throw std::runtime_error("AlgebraicOperator<DGVertex>::add_exit_arc -- trying to add an arc to a vertex not equivalent to either argument.");
  }
}

