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

#ifndef _libint2_src_bin_libint_entity_h_
#define _libint2_src_bin_libint_entity_h_

#include <class_registry.h>
#include <dgvertex.h>

#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

namespace libint2 {

/**
   EntityTypes enumerates the types of objects Entity can represent
*/
namespace EntityTypes {
typedef enum { fp, integer } EntityTypeEnum;

template <unsigned int TypeIndex>
struct EntityType {
  static unsigned int type2int() { return TypeIndex; }
};

typedef EntityType<fp> FP;
typedef EntityType<integer> Int;

static const unsigned int ntypes = 2;
};  // namespace EntityTypes

/// Product of 2 types
template <typename T, typename U>
struct ProductType;
template <>
struct ProductType<int, int> {
  typedef int result;
};
template <>
struct ProductType<int, double> {
  typedef double result;
};
template <>
struct ProductType<double, int> {
  typedef double result;
};
template <>
struct ProductType<double, double> {
  typedef double result;
};
template <>
struct ProductType<EntityTypes::Int, int> {
  typedef EntityTypes::Int result;
};
template <>
struct ProductType<EntityTypes::Int, double> {
  typedef EntityTypes::FP result;
};
template <>
struct ProductType<EntityTypes::FP, int> {
  typedef EntityTypes::FP result;
};
template <>
struct ProductType<EntityTypes::FP, double> {
  typedef EntityTypes::FP result;
};
template <>
struct ProductType<int, EntityTypes::Int> {
  typedef EntityTypes::Int result;
};
template <>
struct ProductType<double, EntityTypes::Int> {
  typedef EntityTypes::FP result;
};
template <>
struct ProductType<int, EntityTypes::FP> {
  typedef EntityTypes::FP result;
};
template <>
struct ProductType<double, EntityTypes::FP> {
  typedef EntityTypes::FP result;
};
template <>
struct ProductType<EntityTypes::Int, EntityTypes::Int> {
  typedef EntityTypes::Int result;
};
template <>
struct ProductType<EntityTypes::Int, EntityTypes::FP> {
  typedef EntityTypes::FP result;
};
template <>
struct ProductType<EntityTypes::FP, EntityTypes::Int> {
  typedef EntityTypes::FP result;
};
template <>
struct ProductType<EntityTypes::FP, EntityTypes::FP> {
  typedef EntityTypes::FP result;
};

/// Converts x to its string representation
template <typename T>
std::string to_string(const T& x) {
  std::ostringstream oss;
  oss << std::scientific
      << std::setprecision(std::numeric_limits<T>::digits10 + 1) << x;
  return oss.str();
}

/**
Entity is a base class for all objects that exist at compile or runtime of the
generated code.
*/
class Entity {
 public:
  virtual ~Entity() {}
  /// Return id string
  const std::string& id() const { return id_; }

 protected:
  Entity(const std::string& id) : id_(id) {}

 private:
  /// Short id label
  std::string id_;
};

/**
   RTimeEntity is an Entity of type T that exists at runtime of the generated
   code (hence has no value known at compile-time)
*/
template <class T>
class RTimeEntity : public Entity, public DGVertex {
 public:
  typedef typename DGVertex::KeyType key_type;

  RTimeEntity(const std::string& id, bool p = true)
      : Entity(id),
        DGVertex(ClassInfo<RTimeEntity>::Instance().id()),
        precomputed_(p) {
    FNVStringHash SH;
    key_ = KeyTypes::cast(SH.hash(id));
#if DEBUG
    std::cout << "Allocated RTimeEntity id = " << this->id() << std::endl;
#endif
  }

  virtual ~RTimeEntity() {
#if DEBUG
    std::cout << "Deallocated RTimeEntity id = " << this->id() << std::endl;
#endif
  }

  /// Implementation of DGVertex::size()
  unsigned int size() const override { return 1; }

  /// Implementation of DGVertex::equiv()
  bool equiv(const std::shared_ptr<DGVertex>& a) const override {
    if (a->typeid_ == typeid_) {
#if USE_INT_KEY_TO_COMPARE
      return key() == a->key() && label() == a->label();
#else
      std::shared_ptr<RTimeEntity> a_cast =
          std::static_pointer_cast<RTimeEntity, DGVertex>(a);
      return id() == a_cast->id();
#endif
    } else
      return false;
  }

  /// Implementation of DGVertex::label()
  const std::string& label() const override { return Entity::id(); }
  /// Implementation of DGVertex::id()
  const std::string& id() const override { return label(); }
  /// Implementation of DGVertex::description()
  std::string description() const override {
    std::ostringstream os;
    os << "RTimeEntity: " << id();
    const std::string descr = os.str();
    return descr;
  }
  /// Implements Hashable::key()
  typename DGVertex::KeyReturnType key() const override { return key_; }

 private:
  /// Implementation of DGVertex::this_precomputed()
  bool this_precomputed() const override { return precomputed_; }

  key_type key_;
  /// RTimeEntity can go either way. Example of a precomputed quartity is a
  /// quantity passed by the user to the library. Example of a quantity that is
  /// not precomputed and thus must be evaluated explicitly by the graph is a
  /// product of precomputed RTimeEntity with a CTimeEntity
  bool precomputed_;
};

/**
   CTimeEntity is an Entity of type T that exists at compile-time of the
   generated code (hence has a value known at compile-time)
*/
template <class T>
class CTimeEntity : public Entity, public DGVertex {
 public:
  CTimeEntity(const T& val)
      : Entity(to_string(val)),
        DGVertex(ClassInfo<CTimeEntity>::Instance().id()),
        value_(val) {
#if DEBUG
    std::cout << "Allocated CTimeEntity id = " << this->id()
              << " value = " << value() << std::endl;
#endif
  }

  virtual ~CTimeEntity() {
#if DEBUG
    std::cout << "Deallocated CTimeEntity id = " << this->id()
              << " value = " << value() << std::endl;
#endif
  }

  /// Implementation of DGVertex::size()
  unsigned int size() const override { return 1; }

  /// Implementation of DGVertex::equiv()
  bool equiv(const std::shared_ptr<DGVertex>& a) const override {
    if (a->typeid_ == typeid_) {
#if USE_INT_KEY_TO_COMPARE
      return key() == a->key();
#else
      std::shared_ptr<CTimeEntity> a_cast =
          std::static_pointer_cast<CTimeEntity, DGVertex>(a);
      return id() == a_cast->id();
#endif
    } else
      return false;
  }

  /// Implementation of DGVertex::label()
  const std::string& label() const override { return Entity::id(); }
  /// Implementation of DGVertex::id()
  const std::string& id() const override { return label(); }
  /// Implementation of DGVertex::description()
  std::string description() const override {
    std::ostringstream os;
    os << "CTimeEntity: " << id();
    const std::string descr = os.str();
    return descr;
  }

  /// returns the value
  typename KeyTraits<T>::ReturnType value() const { return value_; }

  /// Implements Hashable::key()
  typename DGVertex::KeyReturnType key() const override {
    if (std::is_floating_point<T>::value) {
      if (not std::is_same<T, double>::value)
        throw std::runtime_error(
            "CTimeEntity<Real> only supported when Real==double");
      return static_cast<typename DGVertex::KeyReturnType>(
          *reinterpret_cast<const unsigned long*>(&value_));
    } else
      return static_cast<typename DGVertex::KeyReturnType>(value());
  }

 private:
  T value_;

  /// Implementation of DGVertex::this_precomputed()
  bool this_precomputed() const override { return true; }
};

/** Creates product A*B. Exact type depends on type of A -- if A is a
   runtime-entity, then the result is a runtime entity as well. Otherwise the
   result is a compile-time entity.
*/
//  template <typename T>
//    std::shared_ptr<Entity>
//    operator*(const std::shared_ptr<Entity>& A, const std::shared_ptr<
//    CTimeEntity<T> >& B);

/** Creates product A*B.
 */
template <typename T, typename U>
std::shared_ptr<CTimeEntity<typename ProductType<T, U>::result> > operator*(
    const std::shared_ptr<CTimeEntity<T> >& A,
    const std::shared_ptr<CTimeEntity<U> >& B) {
  typedef CTimeEntity<typename ProductType<T, U>::result> prodtype;
  return std::shared_ptr<prodtype>(new prodtype(A->value() * B->value()));
}

/** Creates product A*B.
 */
template <typename T, typename U>
std::shared_ptr<RTimeEntity<typename ProductType<T, U>::result> > operator*(
    const std::shared_ptr<RTimeEntity<T> >& A,
    const std::shared_ptr<CTimeEntity<U> >& B) {
  typedef RTimeEntity<typename ProductType<T, U>::result> prodtype;
  std::ostringstream oss;
  oss << A->id() << "*" << B->id();
  // TODO this should be false, but the logic of DirectedGraph construction
  // depends on this being true
  const bool not_precomputed = true;
  return std::shared_ptr<prodtype>(new prodtype(oss.str(), not_precomputed));
}
// TODO should be possible to enable this, but this creates RTimeEntities that
// should not be precomputed, see the comment above
#if 0
  /** Creates product B*A.
  */
  template <typename T, typename U>
    std::shared_ptr< RTimeEntity< typename ProductType<T,U>::result > >
    operator*(const std::shared_ptr< CTimeEntity<U> >& B, const std::shared_ptr< RTimeEntity<T> >& A)
    {
      return A * B;
    }
#endif

};  // namespace libint2

#endif
