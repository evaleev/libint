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

#include <graph_registry.h>

using namespace libint2;

GraphRegistry::GraphRegistry()
    : accumulate_targets_(false),
      return_targets_(true),
      unroll_threshold_(0),
      uncontract_(false),
      ignore_missing_prereqs_(false),
      do_cse_(false),
      condense_expr_(false),
      stack_name_("inteval->stack"),
      current_timer_(-1) {}

GraphRegistry::~GraphRegistry() {}

GraphRegistry* GraphRegistry::clone() const {
  GraphRegistry* gr = new GraphRegistry;
  *gr = *this;
  return gr;
}

////

InternalGraphRegistry::InternalGraphRegistry()
    : accumulate_targets_directly_(false), size_of_target_accum_(0) {}

InternalGraphRegistry::~InternalGraphRegistry() {}
