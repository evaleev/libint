
#include <string>
#include <rr.h>

#ifndef _libint2_src_bin_libint_entity_h_
#define _libint2_src_bin_libint_entity_h_

namespace libint2 {
  
  /** 
  Entity is a base class for all objects that exist at compile or runtime of the generated code.
  */
  class Entity
  {
    
    public:
    virtual ~Entity() {}
    /// Return id string
    std::string id() const { return id_; }
    
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
        Entity(id),DGVertex()
        {
        }

      ~RTimeEntity()
        {
        }

      /// Implementation of DGVertex::size()
      const unsigned int size() const { return 0; }
    
      /// Implementation of DGVertex::equiv()
      bool equiv(const SafePtr<DGVertex>& a) const
      {
	SafePtr<RTimeEntity> a_cast = dynamic_pointer_cast<RTimeEntity,DGVertex>(a);
	if (a_cast) {
	  return id() == a_cast->id();
	}
	else
	  return false;
      }
      
      /// Implementation of DGVertex::print()
      void print(std::ostream& os) const
      {
	os << "RTimeEntity: " << id();
      }
      
      /// Implementation of DGVertex::precomputed()
      bool precomputed() const
      {
        return true;
      }
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
        Entity(id), DGVertex(), value_(val)
        {
        }

      ~CTimeEntity()
        {
        }

      /// Implementation of DGVertex::size()
      const unsigned int size() const { return 0; }

      /// Implementation of DGVertex::equiv()
      bool equiv(const SafePtr<DGVertex>& a) const
      {
	SafePtr<CTimeEntity> a_cast = dynamic_pointer_cast<CTimeEntity,DGVertex>(a);
	if (a_cast) {
	  return id() == a_cast->id();
	}
	else
	  return false;
      }
      
      /// Implementation of DGVertex::print()
      void print(std::ostream& os) const
      {
	os << "CTimeEntity: " << id() << " value = " << value_;
      }
      
      /// Implementation of DGVertex::precomputed()
      bool precomputed() const
      {
        return true;
      }

    private:
      T value_;

    };

};

#endif

