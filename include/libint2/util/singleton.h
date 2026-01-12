/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint library.
 *
 *  Libint library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_include_libint2_util_singleton_h_
#define _libint2_include_libint2_util_singleton_h_

#include <memory>

namespace libint2 {
namespace detail {

template <typename T>
class managed_singleton {
 public:
  static T* instance() {
    if (not instance_) instance_ = std::unique_ptr<T>(new T);
    return instance_.get();
  }
  static bool instance_exists() { return instance_.get() != nullptr; }
  static void delete_instance() { instance_.reset(); }

 private:
  managed_singleton() = delete;
  static std::unique_ptr<T> instance_;
};
template <typename T>
std::unique_ptr<T> managed_singleton<T>::instance_;

}  // namespace detail
}  // namespace libint2

#endif /* header guard */
