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
#include <graph_registry.h>
#include <smart_ptr.h>
#include <tactic.h>

#ifndef _libint2_src_bin_libint_strategy_h_
#define _libint2_src_bin_libint_strategy_h_

#define USE_HRR_FOR_TiG12 1

namespace libint2 {

class DGVertex;
class DirectedGraph;

/**
Strategy specifies how to apply recurrence relations.
*/
class Strategy {
 public:
  typedef std::shared_ptr<RecurrenceRelation> RR;
  Strategy() {}
  ~Strategy() {}

  /// Returns the optimal recurrence relation for integral
  RR optimal_rr(const std::shared_ptr<DirectedGraph>& graph,
                const std::shared_ptr<DGVertex>& integral,
                const std::shared_ptr<Tactic>& tactic);
};

};  // namespace libint2

#endif
