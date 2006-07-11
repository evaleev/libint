
#include <string>
#include <dgvertex.h>

#ifndef _libint2_src_bin_libint_entity_h_
#define _libint2_src_bin_libint_entity_h_

namespace libint2 {

  /**
     EntityTypes enumerates the types of objects Entity can represent
  */
  namespace EntityTypes {
    typedef enum {fp, integer} EntityTypeEnum;

    template <unsigned int TypeIndex>
    struct EntityType {
      static unsigned int type2int() {
        return TypeIndex;
      }
    };

    typedef EntityType<fp> FP;
    typedef EntityType<integer> Int;

    static const unsigned int ntypes = 2;
  };

  /// Product of 2 types
  template <typename T, typename U>
    struct ProductType;
  template <> struct ProductType<int,int>         { typedef int result; };
  template <> struct ProductType<int,double>      { typedef double result; };
  template <> struct ProductType<double,int>      { typedef double result; };
  template <> struct ProductType<double,double>   { typedef double result; };
  template <> struct ProductType<EntityTypes::Int,int>     { typedef EntityTypes::Int result; };
  template <> struct ProductType<EntityTypes::Int,double>  { typedef EntityTypes::FP result; };
  template <> struct ProductType<EntityTypes::FP,int>      { typedef EntityTypes::FP result; };
  template <> struct ProductType<EntityTypes::FP,double>   { typedef EntityTypes::FP result; };
  template <> struct ProductType<int,EntityTypes::Int>     { typedef EntityTypes::Int result; };
  template <> struct ProductType<double,EntityTypes::Int>  { typedef EntityTypes::FP result; };
  template <> struct ProductType<int,EntityTypes::FP>      { typedef EntityTypes::FP result; };
  template <> struct ProductType<double,EntityTypes::FP>   { typedef EntityTypes::FP result; };
  template <> struct ProductType<EntityTypes::Int,EntityTypes::Int>     { typedef EntityTypes::Int result; };
  template <> struct ProductType<EntityTypes::Int,EntityTypes::FP>      { typedef EntityTypes::FP result; };
  template <> struct ProductType<EntityTypes::FP,EntityTypes::Int>      { typedef EntityTypes::FP result; };
  template <> struct ProductType<EntityTypes::FP,EntityTypes::FP>       { typedef EntityTypes::FP result; };

  /// Converts x to its string representation
  template <typename T>
    std::string to_string(const T& x) {
    std::ostringstream oss;  oss << x;  return oss.str();
  }
  
  /** 
  Entity is a base class for all objects that exist at compile or runtime of the generated code.
  */
  class Entity
  {
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
     RTimeEntity is an Entity of type T that exists at runtime of the generated code (hence
     has no value known at compile-time)
  */
  template <class T>
    class RTimeEntity :
    public Entity,
    public DGVertex
    {
      public:
      RTimeEntity(const std::string& id) :
        Entity(id),DGVertex(ClassInfo<RTimeEntity>::Instance().id()), descr_()
        {
          FNVStringHash SH;
          key_ = SH.hash(id);
#if DEBUG
          std::cout << "Allocated RTimeEntity id = " << this->id() << std::endl;
#endif
        }

      ~RTimeEntity()
        {
#if DEBUG
          std::cout << "Deallocated RTimeEntity id = " << this->id() << std::endl;
#endif
        }

      /// Implementation of DGVertex::size()
      const unsigned int size() const { return 1; }
    
      /// Implementation of DGVertex::equiv()
      bool equiv(const SafePtr<DGVertex>& a) const
      {
	if (a->typeid_ == typeid_) {
#if USE_INT_KEY_TO_COMPARE
          return key() == a->key() && label() == a->label();
#else
          SafePtr<RTimeEntity> a_cast = static_pointer_cast<RTimeEntity,DGVertex>(a);
	  return id() == a_cast->id();
#endif
	}
	else
	  return false;
      }
      
      /// Implementation of DGVertex::label()
      const std::string& label() const
      {
	return Entity::id();
      }
      /// Implementation of DGVertex::id()
      const std::string& id() const
      {
	return label();
      }
      /// Implementation of DGVertex::description()
      const std::string& description() const
      {
        if (descr_.empty()) {
          ostringstream os;
          os << "RTimeEntity: " << id();
          descr_ = os.str();
        }
        return descr_;
      }
      /// Implements Hashable::key()
      typename DGVertex::KeyReturnType key() const {
        return key_;
      }
      
      private:
      /// Implementation of DGVertex::this_precomputed()
      bool this_precomputed() const
      {
        return true;
      }

      mutable std::string descr_;
      typename FNVStringHash::KeyType key_;
    };

  /**
     CTimeEntity is an Entity of type T that exists at compile-time of the generated code (hence
     has a value known at compile-time)
  */
  template <class T>
    class CTimeEntity :
    public Entity,
    public DGVertex
    {
      public:
      CTimeEntity(const T& val) :
        Entity(to_string(val)), DGVertex(ClassInfo<CTimeEntity>::Instance().id()), value_(val), descr_()
        {
#if DEBUG
          std::cout << "Allocated CTimeEntity id = " << this->id() << " value = " << value() << std::endl;
#endif
        }

      ~CTimeEntity()
        {
#if DEBUG
          std::cout << "Deallocated CTimeEntity id = " << this->id() << " value = " << value() << std::endl;
#endif
        }

      /// Implementation of DGVertex::size()
      const unsigned int size() const { return 1; }

      /// Implementation of DGVertex::equiv()
      bool equiv(const SafePtr<DGVertex>& a) const
      {
	if (a->typeid_ == typeid_) {
#if USE_INT_KEY_TO_COMPARE
          return key() == a->key();
#else
          SafePtr<CTimeEntity> a_cast = static_pointer_cast<CTimeEntity,DGVertex>(a);
	  return id() == a_cast->id();
#endif
	}
	else
	  return false;
      }
      
      /// Implementation of DGVertex::label()
      const std::string& label() const
      {
	return Entity::id();
      }
      /// Implementation of DGVertex::id()
      const std::string& id() const
      {
	return label();
      }
      /// Implementation of DGVertex::description()
      const std::string& description() const
      {
        if (descr_.empty()) {
          ostringstream os;
          os << "CTimeEntity: " << id();
          descr_ = os.str();
        }
        return descr_;
      }
      
      /// returns the value
      const T& value() const { return value_; }

      /// Implements Hashable::key()
      typename DGVertex::KeyReturnType key() const {
        return static_cast<typename DGVertex::KeyReturnType>(value());
      }

    private:
      T value_;

      /// Implementation of DGVertex::this_precomputed()
      bool this_precomputed() const
      {
        return true;
      }

      mutable std::string descr_;

    };
    
  /** Creates product A*B. Exact type depends on type of A -- if A is a runtime-entity,
      then the result is a runtime entity as well. Otherwise the result is a compile-time entity.
  */
  template <typename T>
    SafePtr<Entity>
    operator*(const SafePtr<Entity>& A, const SafePtr< CTimeEntity<T> >& B);

  /** Creates product A*B.
  */
  template <typename T, typename U>
    SafePtr< CTimeEntity< typename ProductType<T,U>::result > >
    operator*(const SafePtr< CTimeEntity<T> >& A, const SafePtr< CTimeEntity<U> >& B)
    {
      typedef CTimeEntity< typename ProductType<T,U>::result > prodtype;
      return SafePtr<prodtype>(new prodtype(A->value()*B->value()));
    }

  /** Creates product A*B.
  */
  template <typename T, typename U>
    SafePtr< RTimeEntity< typename ProductType<T,U>::result > >
    operator*(const SafePtr< RTimeEntity<T> >& A, const SafePtr< CTimeEntity<U> >& B)
    {
      typedef RTimeEntity< typename ProductType<T,U>::result > prodtype;
      ostringstream oss;
      oss << A->id() << "*" << B->id();
      return SafePtr<prodtype>(new prodtype(oss.str()));
    }

};

#endif

