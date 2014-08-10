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

#ifndef _libint2_src_bin_libint_util_h_
#define _libint2_src_bin_libint_util_h_

#include <numeric>
#include <string>
#include <stdexcept>
#include <smart_ptr.h>
#include <util_types.h>

namespace libint2 {
  std::string to_string(FunctionPosition pos);
  
  template <class Target, class Source> SafePtr<Target> require_dynamic_cast(const SafePtr<Source>& s) {
    const SafePtr<Target> t = dynamic_pointer_cast<Target,Source>(s);
    if (t == 0)
      throw std::runtime_error("require_dynamic_cast: dynamic case failed");
    return t;
  }
  template <class Target, class Source> const Target* require_dynamic_cast(const Source* s) {
    const Target* t = dynamic_cast<Target*>(s);
    if (t == 0)
      throw std::runtime_error("require_dynamic_cast: dynamic case failed");
    return t;
  }
} // namespace libint2

#endif /* header guard */
