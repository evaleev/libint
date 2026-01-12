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

#ifndef _libint2_src_bin_libint_classregistry_h_
#define _libint2_src_bin_libint_classregistry_h_

#include <key.h>

namespace libint2 {

/** This is a unique registry of classes. */
class ClassRegistry {
 public:
  typedef KeyTypes::ClassID ClassID;
  static ClassRegistry& Instance();
  ClassID next_id() { return nclasses_++; }

 private:
  ClassRegistry();
  static ClassRegistry* registry_;
  ClassID nclasses_;
};

/** Objects of this type provide limited information about the class at runtime.
   Unlike type_info, these objects don't have to be constructed using an object
   of type T.
 */
template <typename T>
class ClassInfo {
 public:
  typedef ClassRegistry::ClassID ClassID;

  static ClassInfo& Instance() {
    if (!info_) info_ = new ClassInfo;
    return *info_;
  }

  ~ClassInfo() {}

  ClassID id() const { return id_; }

 private:
  ClassInfo() : id_(ClassRegistry::Instance().next_id()) {}

  static ClassInfo* info_;
  ClassID id_;
};

template <typename T>
ClassInfo<T>* ClassInfo<T>::info_;

};  // namespace libint2

#endif
