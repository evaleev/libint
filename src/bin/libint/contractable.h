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

#ifndef _libint2_src_bin_libint_contract_h_
#define _libint2_src_bin_libint_contract_h_

namespace libint2 {

/// use this as a base to add to Derived a "contracted()" attribute
template <typename Derived>
class Contractable {
 public:
  Contractable() : value_(default_value_) {}
  Contractable(const Contractable& source) : value_(source.value_) {}
  Contractable& operator=(const Contractable& source) {
    value_ = source.value_;
    return *this;
  }
  bool contracted() const { return value_; }
  void uncontract() { value_ = false; }
  void contract() { value_ = true; }
  static void set_contracted_default_value(bool dv) { default_value_ = dv; }

 private:
  bool value_;
  static bool default_value_;
};
template <typename Derived>
bool Contractable<Derived>::default_value_ = false;

};  // namespace libint2

#endif
