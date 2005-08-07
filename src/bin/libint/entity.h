
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
        }

      ~RTimeEntity()
        {
        }

      /// Implementation of DGVertex::size()
      const unsigned int size() const { return 1; }
    
      /// Implementation of DGVertex::equiv()
      bool equiv(const SafePtr<DGVertex>& a) const
      {
	if (a->typeid_ == typeid_) {
          SafePtr<RTimeEntity> a_cast = static_pointer_cast<RTimeEntity,DGVertex>(a);
	  return id() == a_cast->id();
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
      CTimeEntity(const std::string& id, const T& val) :
        Entity(id), DGVertex(ClassInfo<CTimeEntity>::Instance().id()), value_(val), descr_()
        {
        }

      ~CTimeEntity()
        {
        }

      /// Implementation of DGVertex::size()
      const unsigned int size() const { return 1; }

      /// Implementation of DGVertex::equiv()
      bool equiv(const SafePtr<DGVertex>& a) const
      {
	if (a->typeid_ == typeid_) {
          SafePtr<CTimeEntity> a_cast = static_pointer_cast<CTimeEntity,DGVertex>(a);
	  return id() == a_cast->id();
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

};

#endif

