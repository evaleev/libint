
#include <vector>
#include <integral.h>
#include <iter.h>

#ifndef _libint2_src_bin_libint_traits_h_
#define _libint2_src_bin_libint_traits_h_

using namespace std;

namespace libint2 {

  /**
    Definition of a generic trait is provided in traits_gen.h
  */

  //
  // StdLibintTraits<CGShell>
  //
  template <>
  void
  StdLibintTraits<CGShell>::init_subobj(const CGShell* cgshell, vector<const CGF*>& cgfs)
  {
    unsigned int am = cgshell->qn();
    unsigned int qn[3];
    for(unsigned int i=0; i<=am; i++) {
      qn[0] = am - i;
      for(unsigned int j=0; j<=i; j++) {
        qn[1] = i - j;
        qn[2] = j;

        cgfs.push_back(new CGF(qn));
      }
    }
  }

  template <>
  void
  StdLibintTraits<CGShell>::dealloc_subobj(vector<const CGF*>& subobj)
  {
    int nelem = subobj.size();
    for(int i=0; i<nelem; i++)
      subobj[i]->~CGF();
  }

  template <class Oper, class BFS, class BraSetType, class KetSetType>
    struct StdLibintTraits< GenIntegralSet<Oper,BFS,BraSetType,KetSetType> >
    {
      typedef GenIntegralSet<Oper,BFS,BraSetType,KetSetType> obj_type;
      typedef typename obj_type::iter_type subobj_type;
      static const unsigned int np = Oper::Properties::np;

      /**
        Order subobjects by iterating over BFS in BraSetType and KetSetType.
        Order of iteration:
          iterate over operators
            iterate over particles
              iterate over it's bra function sets, then ket function sets
        
        */
      static void init_subobj(const obj_type* obj, vector<const subobj_type*>& subobj) {
        
        vector< SubIterator* > siters_inord; // subiterators used to iterate over each set (in the above order)
        vector< vector< SubIterator* > > bra_siters; // subiterators for bra basis function sets (outer vector runs over particle index)
        vector< vector< SubIterator* > > ket_siters; // subiterators for ket basis function sets (outer vector runs over particle index)
        bra_siters.resize(np);
        ket_siters.resize(np);

        // Obtain subiterators in order
        SubIteratorBase< Oper > oper_siter(&obj->oper());
        siters_inord.push_back(&oper_siter);
        int num_suboper = oper_siter.num_iter();
        for(int o=0; o<num_suboper; o++) {
          for(int p=0; p<np; p++) {
            const unsigned int nbra = obj->bra().num_members(p);
            bra_siters[p].resize(nbra);
            for(int i=0; i<nbra; i++) {
              SubIterator* iter = obj->bra().member_subiter(p,i);
              siters_inord.push_back(iter);
              bra_siters[p][i] = iter;
            }

            const unsigned int nket = obj->ket().num_members(p);
            ket_siters[p].resize(nket);
            for(int i=0; i<nket; i++) {
              SubIterator* iter = obj->ket().member_subiter(p,i);
              siters_inord.push_back(iter);
              ket_siters[p][i] = iter;
            }
          }
        }

        const unsigned int niters = siters_inord.size();
        for(int it=0; it<niters; it++)
          siters_inord[it]->init();

        // Now iterate over contents of each subiterator
        bool can_iterate = true;
        while (can_iterate) {

          typename Oper::iter_type oper(oper_siter.elem());
          
          // Construct and initialize bra
          typename BraSetType::iter_type bra;
          for(int p=0; p<np; p++) {
            const unsigned int nbra = bra_siters[p].size();
            for(int i=0; i<nbra; i++)
              bra.set_member(bra_siters[p][i]->pelem(), p, i);
          }

          // Construct and initialize ket
          typename KetSetType::iter_type ket;
          for(int p=0; p<np; p++) {
            const unsigned int nket = ket_siters[p].size();
            for(int i=0; i<nket; i++)
              ket.set_member(ket_siters[p][i]->pelem(), p, i);
          }

          // construct this subobj
          subobj_type* curr_subobj = subobj_type::Instance(oper,bra,ket);
          subobj.push_back(curr_subobj);

          // update subiterators to refer to the next element
          for(int it=niters-1; it>=0; it--) {
            SubIterator& siter = *siters_inord[it];
            ++siter;
            // If next element exists -> break out, else -> rewind it and try next subiterator
            if (siter) {
              break;
            }
            else {
              siter.init();
              // If all subiterators have been exhausted then we are done iterating
              if (it == 0)
                can_iterate = false;
            }
          }

        }
      }

      // Nothing is done here because GenIntegralSet objects are Singleton-like and don't need to be destroyed
      static void dealloc_subobj(vector<const subobj_type*>& subobj) {
      }
    };
  

  template <class BFS>
    struct StdLibintTraits< TwoPRep_11_11<BFS> >
    {
      typedef TwoPRep_11_11<BFS> obj_type;
      typedef typename obj_type::iter_type subobj_type;
      typedef SubIteratorBase< typename TwoPRep_11_11<BFS>::parent_type > parent_siter;

      static void init_subobj(const obj_type* obj, vector<const subobj_type*>& subobj) {

        const unsigned int m = obj->m();

        // Iterate over all SubIteratorBase<GenIntegralSet::iter_type>
        parent_siter gis_siter(obj);
        for(gis_siter.init(); gis_siter; ++gis_siter) {
          const typename TwoPRep_11_11<BFS>::parent_type::iter_type* curr_gis_ptr = gis_siter.elem();
          const typename TwoPRep_11_11<BFS>::iter_type* curr_subobj = TwoPRep_11_11<BFS>::iter_type::Instance(curr_gis_ptr->bra(), curr_gis_ptr->ket(), m);
          subobj.push_back(curr_subobj);

        }
      }

      // Nothing is done here because TwoPRep_11_11 objects are Singleton-like and don't need to be destroyed
      static void dealloc_subobj(vector<const TwoPRep_11_11<typename BFS::iter_type>*>& subobj) {
      }
    };
  
};

#endif

