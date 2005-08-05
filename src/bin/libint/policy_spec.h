
#include <vector>
#include <smart_ptr.h>
#include <iter.h>
#include <policy.h>
#include <integral_decl.h>

#ifndef _libint2_src_bin_libint_policyspec_h_
#define _libint2_src_bin_libint_policyspec_h_

namespace libint2 {

  /*
    Definition of a generic StdLibintTDPolicy is provided in policy.h
  */

  /**
  StdLibintTDPolicy<CGShell>::init_subobj initializes CGFs in canonical order.
   The functions in order are produced using the following C++ loop:

   for(int i=0; i<=am; i++) {
     qn[0] = am - i;
     for(int j=0; j<=i; j++) {
       qn[1] = i - j;
       qn[2] = j;
     }
   }

   where am is the angular momentum of the shell and qn[3] are the x, y, and z
   exponents.
   */

  template <>
  void
  StdLibintTDPolicy<CGShell>::init_subobj(const StdLibintTDPolicy<CGShell>::obj_stype& cgshell,
                                          std::vector<StdLibintTDPolicy<CGShell>::subobj_stype>& cgfs);

  template <>
  void
  StdLibintTDPolicy<CGShell>::dealloc_subobj(std::vector<StdLibintTDPolicy<CGShell>::subobj_stype>& subobj);
  /* source is in policy_spec.cc */

  /** StdLibintTDPolicy<GenIntegralSet> describes how integral sets are composed
      of integrals in canonical order.

      Order integrals by iterating over BFS in BraSetType and KetSetType.
        Order of iteration:
          iterate over operators
            iterate over particles
              iterate over it's bra function sets, then ket function sets
   */

  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    struct StdLibintTDPolicy< GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta> >
    {
      typedef GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta> obj_type;
      typedef typename obj_type::iter_type subobj_type;
      static const unsigned int np = Oper::Properties::np;
      /// how these objects are stored
      typedef typename TypeTraits<obj_type>::StorageType obj_stype;
      /// how these subobjects are stored
      typedef typename TypeTraits<subobj_type>::StorageType subobj_stype;

      static void init_subobj(const SafePtr<obj_type>& obj, std::vector< SafePtr<subobj_type> >& subobj) {
        
        std::vector< SubIterator* > siters_inord; // subiterators used to iterate over each set (in the above order)
        std::vector< std::vector< SubIterator* > > bra_siters; // subiterators for bra basis function sets (outer vector runs over particle index)
        std::vector< std::vector< SubIterator* > > ket_siters; // subiterators for ket basis function sets (outer vector runs over particle index)
        bra_siters.resize(np);
        ket_siters.resize(np);

        // Obtain subiterators in order
        SubIteratorBase< Oper > oper_siter(obj->oper());
        siters_inord.push_back(&oper_siter);

        SubIteratorBase< AuxQuanta > aux_siter(obj->aux());
        siters_inord.push_back(&aux_siter);

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

        const unsigned int niters = siters_inord.size();
        for(int it=0; it<niters; it++)
          siters_inord[it]->init();

        // Now iterate over contents of each subiterator
        bool can_iterate = true;
        while (can_iterate) {

          typename Oper::iter_type oper(oper_siter.elem());
          typename AuxQuanta::iter_type aux(aux_siter.elem());
          
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
	  SafePtr<subobj_type> curr_subobj_sptr = subobj_type::Instance(bra,ket,aux,oper);
          subobj.push_back(curr_subobj_sptr);

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
      static void dealloc_subobj(std::vector< SafePtr<subobj_type> >& subobj) {
      }
    };
  

  /** StdLibintTDPolicy<TwoPRep_11_11> should go away soon.
  */

  template <class BFS>
    struct StdLibintTDPolicy< TwoPRep_11_11<BFS> >
    {
      typedef TwoPRep_11_11<BFS> obj_type;
      typedef typename obj_type::iter_type subobj_type;
      typedef SubIteratorBase< typename TwoPRep_11_11<BFS>::parent_type > parent_siter;
      /// how these objects are stored
      typedef typename TypeTraits<obj_type>::StorageType obj_stype;
      /// how these subobjects are stored
      typedef typename TypeTraits<subobj_type>::StorageType subobj_stype;

      static void init_subobj(const SafePtr<obj_type>& obj, std::vector< SafePtr<subobj_type> >& subobj) {

        // Iterate over all SubIteratorBase<GenIntegralSet::iter_type>
        parent_siter gis_siter(obj);
        for(gis_siter.init(); gis_siter; ++gis_siter) {
          const SafePtr<typename TwoPRep_11_11<BFS>::parent_type::iter_type> curr_gis_ptr = gis_siter.elem();
          const SafePtr<subobj_type> curr_subobj =
            subobj_type::Instance(curr_gis_ptr->bra(), curr_gis_ptr->ket(), *curr_gis_ptr->aux().get());
          subobj.push_back(curr_subobj);
        }
      }

      // Nothing is done here because TwoPRep_11_11 objects are Singleton-like and don't need to be destroyed
      static void dealloc_subobj(std::vector< SafePtr< TwoPRep_11_11<typename BFS::iter_type> > >& subobj) {
      }
    };

  /** StdLibintTDPolicy<R12kG12_11_11> should go away soon.
  */

  template <class BFS, int K>
    struct StdLibintTDPolicy< R12kG12_11_11<BFS,K> >
    {
      typedef R12kG12_11_11<BFS,K> obj_type;
      typedef typename obj_type::iter_type subobj_type;
      typedef SubIteratorBase< typename obj_type::parent_type > parent_siter;
      /// how these objects are stored
      typedef typename TypeTraits<obj_type>::StorageType obj_stype;
      /// how these subobjects are stored
      typedef typename TypeTraits<subobj_type>::StorageType subobj_stype;

      static void init_subobj(const SafePtr<obj_type>& obj, std::vector< SafePtr<subobj_type> >& subobj) {

        // Iterate over all SubIteratorBase<GenIntegralSet::iter_type>
        parent_siter gis_siter(obj);
        for(gis_siter.init(); gis_siter; ++gis_siter) {
          const SafePtr<typename obj_type::parent_type::iter_type> curr_gis_ptr = gis_siter.elem();
          const SafePtr<subobj_type> curr_subobj =
            subobj_type::Instance(curr_gis_ptr->bra(), curr_gis_ptr->ket(), *curr_gis_ptr->aux().get());
          subobj.push_back(curr_subobj);
        }
      }

      // Nothing is done here because TwoPRep_11_11 objects are Singleton-like and don't need to be destroyed
      static void dealloc_subobj(std::vector< SafePtr< R12kG12_11_11<typename BFS::iter_type,K> > >& subobj) {
      }
    };

  
  /** StdLibintTDPolicy<TiG12_11_11> should go away soon.
  */

  template <class BFS, int K>
    struct StdLibintTDPolicy< TiG12_11_11<BFS,K> >
    {
      typedef TiG12_11_11<BFS,K> obj_type;
      typedef typename obj_type::iter_type subobj_type;
      typedef SubIteratorBase< typename obj_type::parent_type > parent_siter;
      /// how these objects are stored
      typedef typename TypeTraits<obj_type>::StorageType obj_stype;
      /// how these subobjects are stored
      typedef typename TypeTraits<subobj_type>::StorageType subobj_stype;

      static void init_subobj(const SafePtr<obj_type>& obj, std::vector< SafePtr<subobj_type> >& subobj) {

        // Iterate over all SubIteratorBase<GenIntegralSet::iter_type>
        parent_siter gis_siter(obj);
        for(gis_siter.init(); gis_siter; ++gis_siter) {
          const SafePtr<typename obj_type::parent_type::iter_type> curr_gis_ptr = gis_siter.elem();
          const SafePtr<subobj_type> curr_subobj =
            subobj_type::Instance(curr_gis_ptr->bra(), curr_gis_ptr->ket(), *curr_gis_ptr->aux().get());
          subobj.push_back(curr_subobj);
        }
      }

      // Nothing is done here because TwoPRep_11_11 objects are Singleton-like and don't need to be destroyed
      static void dealloc_subobj(std::vector< SafePtr< TiG12_11_11<typename BFS::iter_type,K> > >& subobj) {
      }
    };
};

#endif

