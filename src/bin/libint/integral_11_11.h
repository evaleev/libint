
#ifndef _libint2_src_bin_libint_integral1111_h_
#define _libint2_src_bin_libint_integral1111_h_

#include <integral.h>

using namespace std;

namespace libint2 {
  
  /**
     Generic integral over a two-body operator with one bfs for each particle in bra and ket.
  */
  template <class BFS, class Oper, class AuxQuanta = EmptySet> class GenIntegralSet_11_11 :
    public GenIntegralSet< Oper, IncableBFSet, typename DefaultTwoPBraket<BFS>::Result, typename DefaultTwoPBraket<BFS>::Result, AuxQuanta >
    {
    public:
      typedef BFS BasisFunctionType;
      typedef Oper OperType;
      typedef typename DefaultTwoPBraket<BFS>::Result BraType;
      typedef typename DefaultTwoPBraket<BFS>::Result KetType;
      typedef AuxQuanta AuxIndexType;
      typedef GenIntegralSet_11_11<BFS,Oper,AuxQuanta> this_type;
      typedef GenIntegralSet< OperType, IncableBFSet, BraType, KetType, AuxIndexType > parent_type;
      
      /// this is a set of these subobjects
      typedef GenIntegralSet_11_11<typename BFS::iter_type, typename Oper::iter_type, typename AuxQuanta::iter_type> iter_type;
      typedef typename parent_type::key_type key_type;
      /// This the type of the object that manages objects of this type as Singletons
      typedef SingletonStack<this_type,key_type> SingletonManagerType;
      /// This class provides comparison operations on pointers
      typedef PtrEquiv<this_type> PtrComp;

      /** This "constructor" takes basis function sets, in Mulliken ordering.
          Returns a pointer to a unique instance, a la Singleton.
          Note that the ordering of arguments is a bit counterintuitive,
          but in fact corresponds to their practical (rather than logical) importance.
      */
      static const SafePtr<this_type> Instance(const BFS& bra0, const BFS& ket0, const BFS& bra1, const BFS& ket1, const AuxIndexType& aux = AuxIndexType(), const OperType& oper = OperType());
      /** This "constructor" uses a wedge of 2 physicists brakets.
      */
      static const SafePtr<this_type>
      Instance(const algebra::Wedge< BraketPair<BFS,PBra>, BraketPair<BFS,PKet> >& braket_wedge,
               const AuxIndexType& aux = AuxIndexType(), const OperType& oper = OperType()) {
        return Instance(braket_wedge.left[0],
                        braket_wedge.right[0],
                        braket_wedge.left[1],
                        braket_wedge.right[1],
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
      /** Returns a pointer to a unique instance, a la Singleton.
          Note that the ordering of arguments is a bit counterintuitive,
          but in fact corresponds to their practical (rather than logical) importance.
      */
      static const SafePtr<this_type> Instance(const BraType& bra, const KetType& ket, const AuxIndexType& aux = AuxIndexType(), const OperType& oper = OperType());
      ~GenIntegralSet_11_11();

      /// Comparison operator
      bool operator==(const this_type&) const;

    private:
      /// This constructor is also private and not implemented since all Integral's are Singletons. Use Instance instead.
      GenIntegralSet_11_11(const OperType& oper, const BraType& bra, const KetType& ket, const AuxIndexType& aux);

      // This is used to manage GenIntegralSet objects as singletons
      static SingletonManagerType singl_manager_;

      /// Implements DGVertex::this_precomputed()
      bool this_precomputed() const;

    };

#if USE_INT_KEY_TO_HASH
  template <class BFS, class Oper, class AuxQuanta>
    typename GenIntegralSet_11_11<BFS,Oper,AuxQuanta>::SingletonManagerType
    GenIntegralSet_11_11<BFS,Oper,AuxQuanta>::singl_manager_(&this_type::key);
#else
#  error "USE_INT_KEY_TO_HASH must be set"
#endif
  
  template <class BFS, class Oper, class AuxQuanta>
    GenIntegralSet_11_11<BFS,Oper,AuxQuanta>::GenIntegralSet_11_11(const OperType& oper, const BraType& bra, const KetType& ket,  const AuxIndexType& aux) :
    parent_type(oper, bra, ket, aux)
    {
      if (bra.num_members(0) != 1)
        throw std::runtime_error("GenIntegralSet_11_11::GenIntegralSet_11_11(bra,ket) -- number of BFSs in bra for particle 0 must be 1");
      if (bra.num_members(1) != 1)
        throw std::runtime_error("GenIntegralSet_11_11::GenIntegralSet_11_11(bra,ket) -- number of BFSs in bra for particle 1 must be 1");
      if (ket.num_members(0) != 1)
        throw std::runtime_error("GenIntegralSet_11_11::GenIntegralSet_11_11(bra,ket) -- number of BFSs in ket for particle 0 must be 1");
      if (ket.num_members(1) != 1)
        throw std::runtime_error("GenIntegralSet_11_11::GenIntegralSet_11_11(bra,ket) -- number of BFSs in ket for particle 1 must be 1");
#if DEBUG
      std::cout << "GenIntegralSet_11_11: constructed " << this->label() << std::endl;
#endif
    }

  template <class BFS, class Oper, class AuxQuanta>
    GenIntegralSet_11_11<BFS,Oper,AuxQuanta>::~GenIntegralSet_11_11()
    {
#if DEBUG
      std::cout << "GenIntegralSet_11_11: destructed " << this->label() << std::endl;
#endif
    }

  template <class BFS, class Oper, class AuxQuanta>
    const SafePtr< GenIntegralSet_11_11<BFS,Oper,AuxQuanta> >
    GenIntegralSet_11_11<BFS,Oper,AuxQuanta>::Instance(const BraType& bra, const KetType& ket,
                                                       const AuxIndexType& aux, const OperType& oper)
    {
      typedef typename SingletonManagerType::value_type map_value_type;
      key_type key = compute_key(oper,bra,ket,aux);
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
    const SafePtr< GenIntegralSet_11_11<BFS,Oper,AuxQuanta> >
    GenIntegralSet_11_11<BFS,Oper,AuxQuanta>::Instance(const BFS& bra0, const BFS& ket0, const BFS& bra1, const BFS& ket1,
                                                       const AuxIndexType& aux, const OperType& oper)
    {
#if USE_BRAKET_H
      typedef BFS BFSRef;
      BFSRef bra0_ref(bra0);
      BFSRef bra1_ref(bra1);
      BFSRef ket0_ref(ket0);
      BFSRef ket1_ref(ket1);
#else
      typedef SafePtr<BFS> BFSRef;
      BFSRef bra0_ref(new BFS(bra0));
      BFSRef bra1_ref(new BFS(bra1));
      BFSRef ket0_ref(new BFS(ket0));
      BFSRef ket1_ref(new BFS(ket1));
#endif
      vector<BFSRef> vbra0;  vbra0.push_back(bra0_ref);
      vector<BFSRef> vbra1;  vbra1.push_back(bra1_ref);
      vector<BFSRef> vket0;  vket0.push_back(ket0_ref);
      vector<BFSRef> vket1;  vket1.push_back(ket1_ref);
      vector< vector<BFSRef> > vvbra;  vvbra.push_back(vbra0);  vvbra.push_back(vbra1);
      vector< vector<BFSRef> > vvket;  vvket.push_back(vket0);  vvket.push_back(vket1);
      BraType bra(vvbra);
      KetType ket(vvket);
      return Instance(bra,ket,aux,oper);
    }
  
  template <class BFS, class Oper, class AuxQuanta>
    bool
    GenIntegralSet_11_11<BFS,Oper,AuxQuanta>::operator==(const this_type& a) const
    {
      return parent_type::PtrComp::equiv(static_cast<const parent_type*>(this),a);
    }

  template <class BFS, class Oper, class AuxQuanta>
    bool
    GenIntegralSet_11_11<BFS,Oper,AuxQuanta>::this_precomputed() const
    {
      return false;
    }

};

#endif

