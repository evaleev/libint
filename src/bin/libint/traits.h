
#include <vector>

#include <policy.h>
#include <integral.h>

#ifndef _libint2_src_bin_libint_traits_h_
#define _libint2_src_bin_libint_traits_h_

using namespace std;

namespace libint2 {

  template < class T>
  struct StdLibintTraits {
    static void init_subobj(const T* cgshell, vector<const typename T::iter_type*>& cgfs);
    static void dealloc_subobj(vector<const typename T::iter_type*>& cgfs);
  };

  template <class T, class P> class SubIteratorBase;
  template < class BFS>
    struct StdLibintTraits< VectorBraket<BFS> >
    {
      static void init_subobj(const VectorBraket<BFS>* obj, vector<const VectorBraket< typename BFS::iter_type>*>& subobj) {

        vector< vector< const SubIteratorBase<BFS,StdLibintPolicy>* > > iters;

        int np = obj->num_part();
        iters.resize(np);
        for(int p=0; p<np; p++) {
          int nf = obj->num_members(p);

          for(int f=0; f<nf; f++) {
            iters[p].push_back(new SubIteratorBase<BFS,StdLibintPolicy>(obj->member(p,f))
                               );
          }
        }
        
      }

      static void dealloc_subobj(vector<const VectorBraket< typename BFS::iter_type>*>& subobj) {
        int nelem = subobj.size();
        for(int i=0; i<nelem; i++)
          subobj[i]->~VectorBraket<typename BFS::iter_type*>();
      }
    };
  
};

#endif
