
#include <vector>
#include <typelist.h>
#include <smart_ptr.h>

#ifndef _libint2_src_bin_libint_policy_h_
#define _libint2_src_bin_libint_policy_h_

using namespace std;


namespace libint2 {

  /**
    ExistsDefaultSubobjAllocator<T,cond> is a helper class that defines the default subobj allocator
    for T. Such allocator is defined only if cond is true.
  */
  template <class T, bool exists>
    struct ExistsDefaultSubobjAllocator;
  
  template <class T>
    struct ExistsDefaultSubobjAllocator<T,true>{
      typedef T obj_type;
      typedef typename obj_type::iter_type subobj_type;

      static void default_init_subobj(const SafePtr<obj_type>& obj, vector< SafePtr<subobj_type> >& subobj)
      {
        subobj.push_back(obj);
      }
      static void default_dealloc_subobj(vector< SafePtr<subobj_type> >& subobj)
      {
      }
    };
      
    
  /** 
      StdLibintTDPolicy<T> is the default type-specific policy.
  */
  template < class T>
    struct StdLibintTDPolicy {
      typedef T obj_type;
      typedef typename obj_type::iter_type subobj_type;

      /// This function allocates subobj of obj (e.g. basis functions contained in a shell)
      static void init_subobj(const SafePtr<obj_type>& obj, vector< SafePtr<subobj_type> >& subobj)
      {
        // If types are not the same then this function should not be used -- user must provide a specialization
        ExistsDefaultSubobjAllocator< T, IsSameType<obj_type,subobj_type>::value >::default_init_subobj(obj,subobj);
      }
      static void dealloc_subobj(vector< SafePtr<subobj_type> >& subobj)
      {
        // If types are not the same then this function should not be used -- user must provide a specialization
        ExistsDefaultSubobjAllocator< T, IsSameType<obj_type,subobj_type>::value >::default_dealloc_subobj(subobj);
      }
    };


  /**
      StdLibintTIPolicy is the default type-independent policy.
   */
  class StdLibintTIPolicy {

    /// This controls whether can unroll an integral set
    static unsigned int max_set_size_to_unroll_;

  public:

    StdLibintTIPolicy() {}

    virtual void set_max_set_size_to_unroll(unsigned int i)
    {
      max_set_size_to_unroll_ = i;
    }

    virtual unsigned int max_set_size_to_unroll() const
    {
      return max_set_size_to_unroll_;
    }
    
  };


  /**
    Policy<T, TDPol, TIPol> defines a policy for type T as a combination of type-independent (TIPol) policies
   and type-specific (TDPol) policies.
   */
  template <class T, class TIPol = StdLibintTIPolicy, template <class> class TDPol = StdLibintTDPolicy>
    class Policy : public TDPol<T>, public TIPol
    {

      /// Returns whether to unroll an IntegralSet
      bool can_unroll_intset(const SafePtr<T>& iset)
      {
        return iset->set_size() <= TIPol::max_set_size_to_unroll();
      }
      
    public:
      typedef typename TDPol<T>::obj_type obj_type;
      typedef typename obj_type::iter_type subobj_type;
      
      static void init_subobj(const SafePtr<obj_type>& obj, vector< SafePtr<subobj_type> >& subobj)
      {
        TDPol<T>::init_subobj(obj,subobj);
      }

      static void dealloc_subobj(const SafePtr<obj_type>& obj, vector< SafePtr<subobj_type> >& subobj)
      {
        TDPol<T>::dealloc_subobj(obj,subobj);
      }
      
  };

  
};

#endif
