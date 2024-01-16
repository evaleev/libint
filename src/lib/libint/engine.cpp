/*
 *  Copyright (C) 2022-2024 Edward F. Valeev
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

#include <libint2/engine.impl.h>

// need to instantiate this just in case user constructs a libint2::Engine using
// default operator params
template __libint2_engine_inline libint2::any
libint2::Engine::enforce_params_type<
   libint2::detail::default_operator_traits::oper_params_type>(
   libint2::Operator oper,
   const libint2::detail::default_operator_traits::oper_params_type& params,
   bool throw_if_wrong_type);
