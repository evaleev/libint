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

#include <smart_ptr.h>
#include <tactic.h>
#include <global_macros.h>
#include <graph_registry.h>

#ifndef _libint2_src_bin_libint_strategy_h_
#define _libint2_src_bin_libint_strategy_h_

#define USE_HRR_FOR_TiG12 1

using namespace std;


namespace libint2 {

  class DGVertex;
  class DirectedGraph;

  /**
  Strategy specifies how to apply recurrence relations.
  */
  class Strategy {

  public:
    typedef SafePtr<RecurrenceRelation> RR;
    Strategy() {}
    ~Strategy() {}

    /// Returns the optimal recurrence relation for integral
    RR optimal_rr(const SafePtr<DirectedGraph>& graph,
                  const SafePtr<DGVertex>& integral,
                  const SafePtr<Tactic>& tactic);
  };

};

#endif
