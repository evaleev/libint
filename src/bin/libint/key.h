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

#ifndef _libint2_src_bin_libint_key_h_
#define _libint2_src_bin_libint_key_h_

#include <gmpxx.h>
#include <libint2/util/intrinsic_types.h>

#include <cstddef>
#include <sstream>

namespace libint2 {

/// Type/Instance combination serves as a key to quickly compare 2 polymorphic
/// Singletons
template <typename T, typename I>
class TypeAndInstance {
 public:
  typedef T Type;
  typedef I Instance;
  TypeAndInstance() : t_(invalid_type_), i_(invalid_instance_) {}
  TypeAndInstance(const Type& t, const Instance& i) : t_(t), i_(i) {}
  TypeAndInstance(const TypeAndInstance& i) : t_(i.t_), i_(i.i_) {}
  const TypeAndInstance& operator=(const TypeAndInstance& i) {
    t_ = i.t_;
    i_ = i.i_;
    return *this;
  }

  const Type& type() const { return t_; }
  const Instance& instance() const { return i_; }

 private:
  Type t_;
  Instance i_;

  static Type invalid_type_;
  static Instance invalid_instance_;
};

template <typename T, typename I>
typename TypeAndInstance<T, I>::Type TypeAndInstance<T, I>::invalid_type_(-1);
template <typename T, typename I>
typename TypeAndInstance<T, I>::Instance
    TypeAndInstance<T, I>::invalid_instance_(-1);

template <typename T, typename I>
bool operator==(const TypeAndInstance<T, I>& a,
                const TypeAndInstance<T, I>& b) {
  return a.type() == b.type() && a.instance() == b.instance();
}

template <typename T, typename I>
bool operator<(const TypeAndInstance<T, I>& a, const TypeAndInstance<T, I>& b) {
  bool result = (a.type() < b.type()) ||
                ((a.type() == b.type()) && (a.instance() < b.instance()));
  return result;
}

///
/// Collection of types used for constructing keys in libint2
///
struct KeyTypes {
  /// distinct classes have unique ClassID's
  typedef unsigned int ClassID;
  /// some classes need to have distinct instances to have unique InstanceID's,
  /// e.g. generalized Singletons
  typedef mpz_class InstanceID;

 private:
  /// mpz_class cannot be constructed via long long and long double, use
  /// std::string instead
  template <typename Target, typename Source>
  static Target string_cast(Source s) {
    std::ostringstream oss;
    oss << s;
    return Target(oss.str());
  }

 public:
  template <typename U>
  inline static InstanceID cast(U i) {
    return InstanceID(i);
  }
};

template <>
inline KeyTypes::InstanceID KeyTypes::cast<long long>(long long i) {
  return string_cast<InstanceID>(i);
}
template <>
inline KeyTypes::InstanceID KeyTypes::cast<unsigned long long>(
    unsigned long long i) {
  return string_cast<InstanceID>(i);
}

/// this composite hashing key works for DGVertex
typedef TypeAndInstance<KeyTypes::ClassID, KeyTypes::InstanceID> DGVertexKey;

};  // namespace libint2

#endif

// Local Variables:
// mode: c++
// End:
