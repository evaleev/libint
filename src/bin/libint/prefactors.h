
#ifndef _libint2_src_bin_libint_prefactors_h_
#define _libint2_src_bin_libint_prefactors_h_

namespace libint2 {

  /** 
      Entity is a base class for all objects that exist at compile or runtime of the generated code.
   */
  class Entity {

  public:
    virtual ~Entity() {}

  protected:
    Entity() {}

  };


  /**
     RTimeEntity is an Entity of type T that exists at runtime of the generated code (hence
     has no value known at compile-time)
  */
  template <class T>
    class RTimeEntity :
    public Entity
    {

    public:
      RTimeEntity() :
        Entity()
        {
        }

      ~RTimeEntity()
        {
        }

    };

  /**
     CTimeEntity is an Entity of type T that exists at compile-time of the generated code (hence
     has a value known at compile-time)
  */
  template <class T>
    class CTimeEntity :
    public Entity
    {

    public:
      CTimeEntity(const T& val) :
        Entity(), value_(val)
        {
        }

      ~CTimeEntity()
        {
        }

    private:
      T value_;

    };

};

#endif
