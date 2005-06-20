

#ifndef _libint2_src_bin_libint_classregistry_h_
#define _libint2_src_bin_libint_classregistry_h_

namespace libint2 {

  /** This is a unique registry of classes. */
  class ClassRegistry {
  public:
    typedef unsigned int ClassID;
    static ClassRegistry& Instance();
    ClassID next_id() { return nclasses_++; }

  private:
    ClassRegistry();
    static ClassRegistry* registry_;
    ClassID nclasses_;
  };

  /** Objects of this type provide limited information about the class at runtime. Unlike type_info,
      these objects don't have to be constructed using an object of type T. 
   */
  template <typename T> class ClassInfo {
  public:
    typedef ClassRegistry::ClassID ClassID;

    static ClassInfo& Instance()
      {
        if (!info_)
          info_ = new ClassInfo;
        return *info_;
      }
    
    ~ClassInfo()
      {
      }

    ClassID id() const { return id_; }
    
  private:
    ClassInfo() :
      id_(ClassRegistry::Instance().next_id())
      {
      }

    static ClassInfo* info_;
    ClassID id_;
  };

  template <typename T>
    ClassInfo<T>*
    ClassInfo<T>::info_;

};

#endif
