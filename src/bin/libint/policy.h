
#include <typelist.h>
#include <smart_ptr.h>

#ifndef _libint2_src_bin_libint_policy_h_
#define _libint2_src_bin_libint_policy_h_

using namespace std;


namespace libint2 {

  /**
  */
  template <class T, bool exists>
    struct ExistsDefaultStdLibintPolicy;
  
  template <class T>
    struct ExistsDefaultStdLibintPolicy<T,true>{
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
      
    
  /** StdLibintPolicy describes assumptions about orderings, etc. as in Libint version 1.
      StdLibintPolicy<T> describes such assumptions about class T. For many parameters
      there is no need to define specializations

  */

  template < class T>
    struct StdLibintPolicy {
      typedef T obj_type;
      typedef typename obj_type::iter_type subobj_type;

      static void init_subobj(const SafePtr<obj_type>& obj, vector< SafePtr<subobj_type> >& subobj)
      {
        // If types are not the same then this function should not be used -- user must provide a specialization
        ExistsDefaultStdLibintPolicy< T, IsSameType<obj_type,subobj_type>::value >::default_init_subobj(obj,subobj);
      }
      static void dealloc_subobj(vector< SafePtr<subobj_type> >& subobj)
      {
        // If types are not the same then this function should not be used -- user must provide a specialization
        ExistsDefaultStdLibintPolicy< T, IsSameType<obj_type,subobj_type>::value >::default_dealloc_subobj(subobj);
      }
    };

};

#endif
