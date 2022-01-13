/*
 *  Copyright (C) 2004-2021 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_bin_libint_integral11_h_
#define _libint2_src_bin_libint_integral11_h_

#include <integral.h>

namespace libint2 {

  /**
     Generic integral over a one-body operator with one bfs for each particle in bra and ket.
  */
  template <class BFS, class Oper, class AuxQuanta = EmptySet> class GenIntegralSet_1_1 :
    public GenIntegralSet< Oper, IncableBFSet, typename DefaultOnePBraket<BFS>::Result, typename DefaultOnePBraket<BFS>::Result, AuxQuanta >
    {
    public:
      typedef BFS BasisFunctionType;
      typedef Oper OperType;
      typedef typename DefaultOnePBraket<BFS>::Result BraType;
      typedef typename DefaultOnePBraket<BFS>::Result KetType;
      typedef AuxQuanta AuxIndexType;
      typedef GenIntegralSet_1_1<BFS,Oper,AuxQuanta> this_type;
      typedef GenIntegralSet< OperType, IncableBFSet, BraType, KetType, AuxIndexType > parent_type;

      /// this is a set of these subobjects
      typedef GenIntegralSet_1_1<typename BFS::iter_type, typename Oper::iter_type, typename AuxQuanta::iter_type> iter_type;
      typedef typename parent_type::key_type key_type;
      /// This the type of the object that manages objects of this type as Singletons
      typedef SingletonStack<this_type,key_type> SingletonManagerType;
      /// This class provides comparison operations on pointers
      typedef PtrEquiv<this_type> PtrComp;

      /** This "constructor" takes basis function sets.
          Returns a pointer to a unique instance, a la Singleton.
          Note that the ordering of arguments is a bit counterintuitive,
          but in fact corresponds to their practical (rather than logical) importance.
      */
      static const SafePtr<this_type> Instance(const BFS& bra0, const BFS& ket0, const AuxIndexType& aux = AuxIndexType(), const OperType& oper = OperType());

#if 0
      /** This "constructor" uses a wedge of 2 physicists brakets.
      */
      static const SafePtr<this_type>
      Instance(const algebra::Wedge< BraketPair<BFS,PBra>, BraketPair<BFS,PKet> >& braket_wedge,
               const AuxIndexType& aux = AuxIndexType(), const OperType& oper = OperType()) {
        return Instance(braket_wedge.left[0],
                        braket_wedge.right[0],
                        aux,
                        oper);
      }
      /** This "constructor" uses a wedge of 2 chemists brakets.
      */
      static const SafePtr<this_type>
      Instance(const algebra::Wedge< BraketPair<BFS,CBra>, BraketPair<BFS,CKet> >& braket_wedge,
               const AuxIndexType& aux = AuxIndexType(), const OperType& oper = OperType()) {
        return Instance(braket_wedge.left[0],
                        braket_wedge.left[1],
                        braket_wedge.right[0],
                        braket_wedge.right[1],
                        aux,
                        oper);
      }
#endif

      /** Returns a pointer to a unique instance, a la Singleton.
          Note that the ordering of arguments is a bit counterintuitive,
          but in fact corresponds to their practical (rather than logical) importance.
      */
      static const SafePtr<this_type> Instance(const BraType& bra, const KetType& ket, const AuxIndexType& aux = AuxIndexType(), const OperType& oper = OperType());
      virtual ~GenIntegralSet_1_1();

      /// Comparison operator
      bool operator==(const this_type&) const;

      /// Reimplements DGVertex::unregister()
      void unregister() const;

      /// Implements GenIntegralSet::auto_unroll()
      bool auto_unroll() const;

   private:
      /// This constructor is also private and not implemented since all Integral's are Singletons. Use Instance instead.
      GenIntegralSet_1_1(const OperType& oper, const BraType& bra, const KetType& ket, const AuxIndexType& aux);

      // This is used to manage GenIntegralSet objects as singletons
      static SingletonManagerType singl_manager_;

      /// Implements DGVertex::this_precomputed()
      bool this_precomputed() const;

    };

#if USE_INT_KEY_TO_HASH
  template <class BFS, class Oper, class AuxQuanta>
    typename GenIntegralSet_1_1<BFS,Oper,AuxQuanta>::SingletonManagerType
    GenIntegralSet_1_1<BFS,Oper,AuxQuanta>::singl_manager_(&this_type::key);
#else
#  error "USE_INT_KEY_TO_HASH must be set"
#endif

  template <class BFS, class Oper, class AuxQuanta>
    GenIntegralSet_1_1<BFS,Oper,AuxQuanta>::GenIntegralSet_1_1(const OperType& oper, const BraType& bra, const KetType& ket,  const AuxIndexType& aux) :
    parent_type(oper, bra, ket, aux)
    {
      if (bra.num_members(0) != 1)
        throw std::runtime_error("GenIntegralSet_1_1::GenIntegralSet_1_1(bra,ket) -- number of BFSs in bra for particle 0 must be 1");
      if (ket.num_members(0) != 1)
        throw std::runtime_error("GenIntegralSet_1_1::GenIntegralSet_1_1(bra,ket) -- number of BFSs in ket for particle 0 must be 1");
#if DEBUG
      std::cout << "GenIntegralSet_1_1: constructed " << this->label() << std::endl;
#endif
    }

  template <class BFS, class Oper, class AuxQuanta>
    GenIntegralSet_1_1<BFS,Oper,AuxQuanta>::~GenIntegralSet_1_1()
    {
#if DEBUG
      std::cout << "GenIntegralSet_1_1: destructed " << this->label() << std::endl;
#endif
    }

  template <class BFS, class Oper, class AuxQuanta>
    const SafePtr< GenIntegralSet_1_1<BFS,Oper,AuxQuanta> >
    GenIntegralSet_1_1<BFS,Oper,AuxQuanta>::Instance(const BraType& bra, const KetType& ket,
                                                       const AuxIndexType& aux, const OperType& oper)
    {
      typedef typename SingletonManagerType::value_type map_value_type;
      key_type key = parent_type::compute_key(oper,bra,ket,aux);
      const map_value_type& val = singl_manager_.find(key);
      if (!val.second) {
        SafePtr<this_type> this_int(new this_type(oper,bra,ket,aux));
        // Use singl_manager_ to make sure this is a new object of this type
        const typename SingletonManagerType::value_type& val = singl_manager_.find(this_int);
        val.second->instid_ = val.first;
        return val.second;
      }
      return val.second;
    }

  template <class BFS, class Oper, class AuxQuanta>
    const SafePtr< GenIntegralSet_1_1<BFS,Oper,AuxQuanta> >
    GenIntegralSet_1_1<BFS,Oper,AuxQuanta>::Instance(const BFS& bra0, const BFS& ket0,
                                                       const AuxIndexType& aux, const OperType& oper)
    {
#if USE_BRAKET_H
      typedef BFS BFSRef;
      BFSRef bra0_ref(bra0);
      BFSRef ket0_ref(ket0);
#else
      typedef SafePtr<BFS> BFSRef;
      BFSRef bra0_ref(new BFS(bra0));
      BFSRef ket0_ref(new BFS(ket0));
#endif
      std::vector<BFSRef> vbra0;  vbra0.push_back(bra0_ref);
      std::vector<BFSRef> vket0;  vket0.push_back(ket0_ref);
      std::vector< std::vector<BFSRef> > vvbra;  vvbra.push_back(vbra0);
      std::vector< std::vector<BFSRef> > vvket;  vvket.push_back(vket0);
      BraType bra(vvbra);
      KetType ket(vvket);
      return Instance(bra,ket,aux,oper);
    }

  template <class BFS, class Oper, class AuxQuanta>
    bool
    GenIntegralSet_1_1<BFS,Oper,AuxQuanta>::operator==(const this_type& a) const
    {
      return parent_type::PtrComp::equiv(static_cast<const parent_type*>(this),a);
    }

  template <class BFS, class Oper, class AuxQuanta>
    void
    GenIntegralSet_1_1<BFS,Oper,AuxQuanta>::unregister() const
    {
      SafePtr<parent_type> this_parent_ptr = const_pointer_cast<parent_type,const parent_type>(EnableSafePtrFromThis<parent_type>::SafePtr_from_this());
      SafePtr<this_type> this_ptr = static_pointer_cast<this_type>(this_parent_ptr);
      singl_manager_.remove(this_ptr);
    }

  // this_precomputed() and auto_unroll() will be specialized, the nonspecialized version is in integral_11_11.impl.h
};

#endif

