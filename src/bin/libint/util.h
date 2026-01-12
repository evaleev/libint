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

#ifndef _libint2_src_bin_libint_util_h_
#define _libint2_src_bin_libint_util_h_

#include <cxxabi.h>
#include <smart_ptr.h>
#include <util_types.h>

#include <numeric>
#include <stdexcept>
#include <string>

namespace libint2 {
std::string to_string(FunctionPosition pos);

template <class Target, class Source>
std::shared_ptr<Target> require_dynamic_cast(const std::shared_ptr<Source>& s) {
  const std::shared_ptr<Target> t =
      std::dynamic_pointer_cast<Target, Source>(s);
  if (t == 0)
    throw std::runtime_error("require_dynamic_cast: dynamic case failed");
  return t;
}
template <class Target, class Source>
const Target* require_dynamic_cast(const Source* s) {
  const Target* t = dynamic_cast<Target*>(s);
  if (t == 0)
    throw std::runtime_error("require_dynamic_cast: dynamic case failed");
  return t;
}

/// @return (demangled) class name
template <typename T>
std::string class_name(T* ptr = nullptr) {
  int status = 1;
  std::unique_ptr<char, void (*)(void*)> result{
      abi::__cxa_demangle(
          ptr == nullptr ? typeid(T).name() : typeid(ptr).name(), NULL, NULL,
          &status),
      std::free};
  return status == 0 ? result.get() : "unknown";
}

}  // namespace libint2

#endif /* header guard */
