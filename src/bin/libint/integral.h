
#include <smart_ptr.h>
#include <dgvertex.h>
#include <oper.h>
#include <typelist.h>
#include <iter.h>
#include <policy_spec.h>
#include <quanta.h>
#include <equiv.h>
#include <singl_stack.h>
#include <global_macros.h>

#if USE_BRAKET_H
#  include <braket.h>
#endif

#ifndef _libint2_src_bin_libint_integral_h_
#define _libint2_src_bin_libint_integral_h_

namespace libint2 {

  using namespace std;
  
  /**
     O is an operator
     TList is a typelist TList1 that consists 2 second-level typelists: TList2 for bra and TList3 for ket
     TList2 and TList3 consist of O::np types which specify basis functions for each particle in bra/ket
  */
  
  // Tags for first and second typelists
  typedef Int2Type<1> Tag1;
  typedef Int2Type<2> Tag2;
  typedef Int2Type<2> Tag3;

  template <class O, class BraList, class KetList>
    class GenIntegral
    {
    public:
      GenIntegral(BraList bra, KetList ket) :
        bra_(bra), ket_(ket) {}
      virtual ~GenIntegral() {}

    private:
      static O O_;
      BraList bra_;
      KetList ket_;

    };


  /**
     This is an abstract base for sets of all types of integrals. Functions can be of any type
     derived from BasisFunctionSet.
  */
  template <class BasisFunctionSet> class IntegralSet {

  public:
    virtual ~IntegralSet() {};

#if USE_BRAKET_H
    /// Obtain pointers to ith BasisFunctionSet for particle p in bra
    virtual const BasisFunctionSet& bra(unsigned int p, unsigned int i) const =0;
    /// Obtain pointers to ith BasisFunctionSet for particle p in ket
    virtual const BasisFunctionSet& ket(unsigned int p, unsigned int i) const =0;
    /// Obtain pointers to ith BasisFunctionSet for particle p in bra
    virtual BasisFunctionSet& bra(unsigned int p, unsigned int i) =0;
    /// Obtain pointers to ith BasisFunctionSet for particle p in ket
    virtual BasisFunctionSet& ket(unsigned int p, unsigned int i) =0;
#else
    /// Obtain pointers to ith BasisFunctionSet for particle p in bra
    virtual const SafePtr<BasisFunctionSet> bra(unsigned int p, unsigned int i) const =0;
    /// Obtain pointers to ith BasisFunctionSet for particle p in ket
    virtual const SafePtr<BasisFunctionSet> ket(unsigned int p, unsigned int i) const =0;
#endif
  };

  /**
     GenIntegralSet is a set of integrals over functions derived from BFS.

     Oper is an operator or a set of operators. Oper must be derived from OperSet.
     BraSetType and KetSetType are types describing sets of functions.
     An example of a class that can be used as BraSetType and KetSetType
     is VectorBraket.
     AuxQuanta describes auxiliary quantum numbers. AuxQuanta should be derived from QuantumSet.
  */
  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta = NullQuantumSet>
    class GenIntegralSet :
    public IntegralSet<BFS>, public DGVertex,
    public EnableSafePtrFromThis< GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta> >
    {
      public:
      /// GenIntegralSet is a set of these subobjects
      typedef GenIntegralSet this_type;
      /// GenIntegralSet is a set of these subobjects
      typedef GenIntegralSet<typename Oper::iter_type, BFS, typename BraSetType::iter_type, typename KetSetType::iter_type, typename AuxQuanta::iter_type> iter_type;
      /// GenIntegralSet is derived from IntegralSet
      typedef IntegralSet<BFS> parent_type;
      /// This type provides comparison operations on pointers to GenIntegralSet
      typedef PtrEquiv<GenIntegralSet> PtrComp;
#if  USE_INT_KEY_TO_HASH
      typedef LIBINT2_UINT_LEAST64 key_type;
      /// This the type of the object that manages GenIntegralSet's as Singletons
      typedef SingletonStack<GenIntegralSet,key_type> SingletonManagerType;
#else
      /// This the type of the object that manages GenIntegralSet's as Singletons
      typedef SingletonStack<GenIntegralSet,std::string> SingletonManagerType;
#endif
      /// This is the type of the operator
      typedef Oper OperType;

      /** No constructors are public since this is a singleton-like quantity.
          Instead, access is provided through Instance().
          
          Derived classes do not have to be Singletons -- hence protected constructors are provided.
      */

      /// GenIntegralSet objects can be handled through the base pointer -- hence the destructor is virtual
      virtual ~GenIntegralSet();

      /// Returns a pointer to a unique instance, a la Singleton
      static const SafePtr<GenIntegralSet> Instance(const BraSetType& bra, const KetSetType& ket, const AuxQuanta& aux, const Oper& oper=Oper());
      
      /// Comparison operator
      virtual bool operator==(const GenIntegralSet&) const;
      /// Specialization of DGVertex::equiv()
      bool equiv(const SafePtr<DGVertex>& v) const
      {
        return PtrComp::equiv(this,v);
      }
      /// Specialization of DGVertex::size()
      virtual const unsigned int size() const;
      /// Specialization of DGVertex::label()
      virtual const std::string& label() const;
      /// Specialization of DGVertex::id()
      virtual const std::string& id() const;
      /// Specialization of DGVertex::description()
      virtual const std::string& description() const;

      /// Implementation of IntegralSet::bra() const
      typename BraSetType::bfs_cref bra(unsigned int p, unsigned int i) const;
      /// Implementation of IntegralSet::ket() const
      typename KetSetType::bfs_cref ket(unsigned int p, unsigned int i) const;
      /// Implementation of IntegralSet::bra()
      typename BraSetType::bfs_ref bra(unsigned int p, unsigned int i);
      /// Implementation of IntegralSet::ket()
      typename KetSetType::bfs_ref ket(unsigned int p, unsigned int i);
      
      typedef BraSetType BraType;
      typedef KetSetType KetType;
      typedef Oper OperatorType;
      typedef AuxQuanta AuxQuantaType;

      /// Obtain the operator
      const SafePtr<Oper> oper() const;
      /// Obtain const ref to bra
      const BraType& bra() const;
      /// Obtain const ref to bra
      const KetType& ket() const;
      /// Obtain the auxiliary quanta
      const SafePtr<AuxQuanta> aux() const;

      /// Implements Hashable::key()
      LIBINT2_UINT_LEAST64 key() const {
        return key_;
      }
      
      /// Reimplements DGVertex::unregister()
      void unregister() const;

      protected:
      // Basic Integral constructor. It is protected so that derived classes don't have to behave like singletons
      GenIntegralSet(const Oper& oper, const BraSetType& bra, const KetSetType& ket, const AuxQuanta& aux);

      BraSetType bra_;
      KetSetType ket_;

      /// set size to sz
      void set_size(unsigned int sz);
      /// Specialization of DGVertex::this_precomputed()
      virtual bool this_precomputed() const { return false; }

      private:
      //
      // All integrals are Singletons by nature, therefore they must be treated as such
      // 1) No public constructors are provided
      // 2) protected members are provided to implement Singleton-type functionality
      //
      GenIntegralSet();
      GenIntegralSet(const GenIntegralSet&);
      // Copy is not permitted
      GenIntegralSet& operator=(const GenIntegralSet& source);

      // This is used to manage GenIntegralSet objects as singletons
      static SingletonManagerType singl_manager_;

      // The operator needs to be a real object rather than real type to be able to construct a SubIterator, etc.
      SafePtr<Oper> O_;
      // Same for AuxQuanta
      SafePtr<AuxQuanta> aux_;

      // size of the integral set
      mutable unsigned int size_;

      // unique label -- no 2 GenIntegralSet instances can have the same label!
      mutable std::string label_;
      // generates label_
      std::string generate_label() const;
      // description
      mutable std::string descr_;
      // key
      mutable key_type key_;

      /// computes and caches key
      void compute_key() const {
        LIBINT2_UINT_LEAST64 key;
        key = ( (O_->key()*bra_.max_key() + bra_.key() ) * ket_.max_key() +
                ket_.key() ) * aux_->max_key() + aux_->key();
        key_ = key;
      }
    };

#if USE_INT_KEY_TO_HASH
  // I use key() to hash GenIntegralSet. Therefore compute_key() must return unique keys!
  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    typename GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::SingletonManagerType
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::singl_manager_(&GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::key);
#else
  // I use label() to hash GenIntegralSet. Therefore labels must be unique!
  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    typename GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::SingletonManagerType
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::singl_manager_(&GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::label);
#endif  

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::GenIntegralSet(const Op& oper, const BraSetType& bra, const KetSetType& ket, const AuxQuanta& aux) :
    DGVertex(ClassInfo<GenIntegralSet>::Instance().id()), O_(SafePtr<Op>(new Op(oper))), bra_(bra), ket_(ket), aux_(SafePtr<AuxQuanta>(new AuxQuanta(aux))),
    size_(0), label_(), descr_()
    {
      if (Op::Properties::np != bra.num_part())
        throw std::runtime_error("GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::GenIntegralSet(bra,ket) -- number of particles in bra doesn't match that in the operator");
      if (Op::Properties::np != ket.num_part())
        throw std::runtime_error("GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::GenIntegralSet(bra,ket) -- number of particles in ket doesn't match that in the operator");
      compute_key();
      std::cout << "Constructed " << label() << std::endl;
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::~GenIntegralSet()
    {
      std::cout << "Destructed " << label() << std::endl;
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const SafePtr< GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta> >
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::Instance(const BraSetType& bra, const KetSetType& ket, const AuxQuanta& aux, const Op& oper)
    {
      SafePtr<this_type> this_int(new this_type(oper,bra,ket,aux));
      // Use singl_manager_ to make sure this is a new object of this type
      const typename SingletonManagerType::value_type& val = singl_manager_.find(this_int);
      val.second->instid_ = val.first;
      this_int.reset();
      return val.second;
    }
  
  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    typename BraSetType::bfs_cref
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::bra(unsigned int p, unsigned int i) const
    {
      return bra_.member(p,i);
    }
    
  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    typename KetSetType::bfs_cref
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::ket(unsigned int p, unsigned int i) const
    {
      return ket_.member(p,i);
    }
    
  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    typename BraSetType::bfs_ref
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::bra(unsigned int p, unsigned int i)
    {
      return bra_.member(p,i);
    }
    
  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    typename KetSetType::bfs_ref
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::ket(unsigned int p, unsigned int i)
    {
      return ket_.member(p,i);
    }
    
  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const typename GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::BraType&
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::bra() const
    {
      return bra_;
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const typename GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::KetType&
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::ket() const
    {
      return ket_;
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const SafePtr<Op>
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::oper() const
    {
      return O_;
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const SafePtr<AuxQuanta>
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::aux() const
    {
      return aux_;
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    bool
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::operator==(const this_type& a) const
    {
#if USE_INT_KEY_TO_COMPARE
      return key() == a.key();
#else
      bool aux_equiv = PtrEquiv<AuxQuanta>::equiv(aux_,a.aux_);
      if (!aux_equiv) return false;
      bool oper_equiv = PtrEquiv<Op>::equiv(O_,a.O_);
      bool bra_equiv = PtrEquiv<BraSetType>::equiv(bra_,a.bra_);
      if (!bra_equiv) return false;
      bool ket_equiv = PtrEquiv<KetSetType>::equiv(ket_,a.ket_);
      if (!ket_equiv) return false;
      return true;
#endif
    }
  
  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    void
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::unregister() const
    {
      singl_manager_.remove(const_pointer_cast<this_type,const this_type>(EnableSafePtrFromThis<this_type>::SafePtr_from_this()));
    }
  
  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const unsigned int
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::size() const
    {
      if (size_ > 0)
        return size_;

#if COMPUTE_SIZE_DIRECTLY
    // WARNING: will not work if aux_ includes "sets" of numbers, i.e. when a quantum number can be a set of numbers
    // but I don't at the moment see why this would be needed
    size_ = bra_.size() * ket_.size() * O_->num_oper();
#else
      // compute size
      SafePtr<this_type> this_ptr = const_pointer_cast<this_type,const this_type>(EnableSafePtrFromThis<GenIntegralSet>::SafePtr_from_this());
      SafePtr< SubIteratorBase<this_type> > siter(new SubIteratorBase<this_type>(this_ptr));
      size_ = siter->num_iter();
      if (size_ == 0)
        throw std::runtime_error("GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::size() -- size is 0");
#endif

      return size_;
    }
    
  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    void
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::set_size(unsigned int sz)
    {
      size_ = sz;
    }
    
  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const std::string&
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::label() const
    {
      if (label_.empty())
        label_ = generate_label();
      return label_;
    }
  
  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    std::string
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::generate_label() const
    {
      ostringstream os;
      os << "< ";
      for(int p=0; p<Op::Properties::np; p++) {
        unsigned int nbra = bra_.num_members(p);
        for(unsigned int i=0; i<nbra; i++)
#if USE_BRAKET_H
          os << bra_.member(p,i).label() << "(" << p << ") ";
#else
          os << bra_.member(p,i)->label() << "(" << p << ") ";
#endif
      }
      os << "| " << O_->label() << " | ";
      for(int p=0; p<Op::Properties::np; p++) {
        unsigned int nket = ket_.num_members(p);
        for(unsigned int i=0; i<nket; i++)
#if USE_BRAKET_H
          os << ket_.member(p,i).label() << "(" << p << ") ";
#else
          os << ket_.member(p,i)->label() << "(" << p << ") ";
#endif
      }
      os << "> ^ { " << aux_->label() << " }";
      return os.str();
    }
  
  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const std::string&
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::id() const
    {
      return label();
    }
  
  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const std::string&
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::description() const
    {
      if (descr_.empty()) {
        ostringstream os;
        os << " GenIntegralSet: " << label();
        descr_ = os.str();
      }
      return descr_;
    }

  /** TwoPRep_11_11_base is the base for all 2-body repulsion integrals with one basis function
    for each particle in bra and ket
    */
  class TwoPRep_11_11_base {
  };

  /// This is the implementation of the Braket concept used by TwoPrep_11_11
  // really need to have typedef template!
  template <typename BFS>
    struct DefaultTwoPBraket {
      /// This defines which Braket implementation to use
      //typedef VectorBraket<BFS> Result;
      typedef ArrayBraket<BFS,2> Result;
    };
  
  /// This is the implementation of the QuantumNumbers concept used by TwoPrep_11_11
  // really need to have typedef template!
  template <typename T, unsigned int N>
    struct DefaultQuantumNumbers {
      /// This defines which QuantumNumbers implementation to use
      //typedef QuantumNumbers<T,N> Result;
      typedef QuantumNumbersA<T,N> Result;
    };
  /**
     mType is the type that describes the auxiliary index of standard 2-body repulsion integrals
  */
  typedef DefaultQuantumNumbers<unsigned int,1>::Result mType;
  /**
     EmptySet is the type that describes null set of auxiliary indices
  */
  typedef DefaultQuantumNumbers<int,0>::Result EmptySet;
  
  /**
     Most basic type -- TwoPRep_11_11 --
     has one bfs for each particle in bra and ket.
     Note that GenIntegralSet is initialized with an abstract type libint2::BFSet,
     from which BFS derives.
  */
  template <class BFS> class TwoPRep_11_11 :
    public GenIntegralSet<TwoERep, IncableBFSet, typename DefaultTwoPBraket<BFS>::Result, typename DefaultTwoPBraket<BFS>::Result, mType >,
    //public EnableSafePtrFromThis< TwoPRep_11_11<BFS> >,
    public TwoPRep_11_11_base
    {
    public:
      typedef TwoERep OperType;
      typedef typename DefaultTwoPBraket<BFS>::Result BraType;
      typedef typename DefaultTwoPBraket<BFS>::Result KetType;
      typedef mType AuxIndexType;
      typedef TwoPRep_11_11 this_type;
      /// TwoPRep_11_11 is a set of these subobjects
      typedef TwoPRep_11_11<typename BFS::iter_type> iter_type;
      /// This is the immediate parent
      typedef GenIntegralSet<TwoERep, IncableBFSet, BraType, KetType, AuxIndexType > parent_type;
      /// This class provides comparison operations on pointers
      typedef PtrEquiv<this_type> PtrComp;

#if USE_INT_KEY_TO_HASH
      typedef typename parent_type::key_type key_type;
      /// This the type of the object that manages GenIntegralSet's as Singletons
      typedef SingletonStack<TwoPRep_11_11,key_type> SingletonManagerType;
#else
      /// This the type of the object that manages GenIntegralSet's as Singletons
      typedef SingletonStack<TwoPRep_11_11,std::string> SingletonManagerType;
#endif
      
      // Destructor only to track object's death
      ~TwoPRep_11_11();
      
      /* This "constructor" takes basis function sets, in Mulliken ordering.
         Returns a pointer to a unique instance, a la Singleton
      */
      static const SafePtr<TwoPRep_11_11> Instance(const BFS& bra0, const BFS& ket0, const BFS& bra1, const BFS& ket1, unsigned int m);
      /// Returns a pointer to a unique instance, a la Singleton
      static const SafePtr<TwoPRep_11_11> Instance(const BraType& bra, const KetType& ket, const AuxIndexType& aux);
      
      unsigned int m() const { return parent_type::aux()->elem(0); };

      /// Comparison operator
      bool operator==(const this_type&) const;
#if OVERLOAD_GENINTEGRALSET_LABEL
      /// Specialization of GenIntegralSet::label()
      const std::string& label() const;
#endif
      /// Reimplements DGVertex::unregister()
      void unregister() const;

    private:
      // This constructor is also private and not implemented since all Integral's are Singletons. Use Instance instead.
      TwoPRep_11_11(const BraType& bra, const KetType& ket, const AuxIndexType& aux);

      // This is used to manage GenIntegralSet objects as singletons
      static SingletonManagerType singl_manager_;

      /// Implements DGVertex::this_precomputed()
      bool this_precomputed() const;
#if OVERLOAD_GENINTEGRALSET_LABEL
      mutable std::string label_;
#endif
    };

#if USE_INT_KEY_TO_HASH
  template <class BFS>
    typename TwoPRep_11_11<BFS>::SingletonManagerType
    TwoPRep_11_11<BFS>::singl_manager_(&TwoPRep_11_11<BFS>::key);
#else
  // I use label() to hash TwoPRep_11_11. Therefore labels must be unique!
  template <class BFS>
    typename TwoPRep_11_11<BFS>::SingletonManagerType
    TwoPRep_11_11<BFS>::singl_manager_(&TwoPRep_11_11<BFS>::label);
#endif
  
  template <class BFS>
    TwoPRep_11_11<BFS>::TwoPRep_11_11(const BraType& bra, const KetType& ket,  const AuxIndexType& aux) :
    GenIntegralSet<TwoERep, IncableBFSet, BraType, KetType, AuxIndexType>(TwoERep(), bra, ket, aux)
    {
      if (bra.num_members(0) != 1)
        throw std::runtime_error("TwoPRep_11_11<BFS>::TwoPRep_11_11(bra,ket) -- number of BFSs in bra for particle 0 must be 1");
      if (bra.num_members(1) != 1)
        throw std::runtime_error("TwoPRep_11_11<BFS>::TwoPRep_11_11(bra,ket) -- number of BFSs in bra for particle 1 must be 1");
      if (ket.num_members(0) != 1)
        throw std::runtime_error("TwoPRep_11_11<BFS>::TwoPRep_11_11(bra,ket) -- number of BFSs in ket for particle 0 must be 1");
      if (ket.num_members(1) != 1)
        throw std::runtime_error("TwoPRep_11_11<BFS>::TwoPRep_11_11(bra,ket) -- number of BFSs in ket for particle 1 must be 1");
      std::cout << "Constructed TwoPRep_11_11 " << this->label() << std::endl;
    }
  
  template <class BFS>
    TwoPRep_11_11<BFS>::~TwoPRep_11_11() {
      std::cout << "Destructed TwoPRep_11_11 " << this->label() << std::endl;
    }
  
  template <class BFS>
    const SafePtr< TwoPRep_11_11<BFS> >
    TwoPRep_11_11<BFS>::Instance(const BraType& bra, const KetType& ket, const AuxIndexType& aux)
    {
      SafePtr<TwoPRep_11_11> this_int(new TwoPRep_11_11<BFS>(bra,ket,aux));
      // Use singl_manager_ to make sure this is a new object of this type
      const typename SingletonManagerType::value_type& val = singl_manager_.find(this_int);
      val.second->instid_ = val.first;
      this_int.reset();
      return val.second;
    }

  template <class BFS>
    const SafePtr< TwoPRep_11_11<BFS> >
    TwoPRep_11_11<BFS>::Instance(const BFS& bra0, const BFS& ket0, const BFS& bra1, const BFS& ket1, unsigned int m)
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
      AuxIndexType aux(vector<unsigned int>(1,m));
      return Instance(bra,ket,aux);
    }
  
  template <class BFS>
    bool
    TwoPRep_11_11<BFS>::operator==(const TwoPRep_11_11<BFS>& a) const
    {
      return parent_type::PtrComp::equiv(static_cast<const parent_type*>(this),a);
    }
  
  template <class BFS>
    void
    TwoPRep_11_11<BFS>::unregister() const
    {
      singl_manager_.remove(
        const_pointer_cast<this_type,const this_type>(
          static_pointer_cast<const this_type, const parent_type>(
            EnableSafePtrFromThis<parent_type>::SafePtr_from_this()
          )
        )
      );
      parent_type::unregister();
    }
  
#if OVERLOAD_GENINTEGRALSET_LABEL
  template <class BFS>
    const std::string&
    TwoPRep_11_11<BFS>::label() const
    {
      if (label_.empty()) {
	ostringstream os;
	os << "(" << parent_type::bra_.member(0,0)->label() << " " << parent_type::ket_.member(0,0)->label()
	   << " | 1/r_{12} | " << parent_type::bra_.member(1,0)->label() << " " << parent_type::ket_.member(1,0)->label() << ")^{" << m() <<"}";
	label_ = os.str();
      }
      return label_;
    };
#endif

  /// TwoPRep_11_11_sq is a shell quartet of ERIs
  typedef TwoPRep_11_11<CGShell> TwoPRep_11_11_sq;

  /// TwoPRep_11_11_int is a single ERIs
  typedef TwoPRep_11_11<CGF> TwoPRep_11_11_int;

  /**
     TypelistBraket is a typelist-based type to describe a bra or a ket
     of an integral. BFS is a base type for all members of a typelist,
     and TList is the typelist itself.
     
     Unfortunately, it seems that the current standard specification does
     not allow for this to work: templates are only processed once, not
     recursively...
  */
  template <class BFS, class TList> class TypelistBraket {
  public:
    typedef TList BFSMatrix;

    TypelistBraket(const TList&);
    TypelistBraket(const TypelistBraket&);
    ~TypelistBraket() throw();

    bool equiv(const TypelistBraket&) const;
    /// Returns pointer to the i-th function for particle p
    const BFS* member(unsigned int p, unsigned int i) const;
    /// Returns the number of BFS for particle p
    const unsigned int num_members(unsigned int p) const;
    /// Returns the number of particles
    const unsigned int num_part() const;

  private:

    BFSMatrix bfs_;

  };

  /*template <class BFS, class TList>
    TypelistBraket<BFS,TList>::TypelistBraket(const TList& bfs) :
    bfs_(bfs)
  {
  }

  template <class BFS, class TList>
    TypelistBraket<BFS,TList>::TypelistBraket(const TypelistBraket<BFS,TList>& a) :
    bfs_(a.bfs_)
  {
  }

  template <class BFS, class TList>
    TypelistBraket<BFS,TList>::~TypelistBraket() throw()
  {
  }

  template <class BFS, class TList>
    const BFS*
    TypelistBraket<BFS,TList>::member(unsigned int p, unsigned int i) const
  {
    return &bfs_.at(p).at(i);
  }
  
  template <class BFS, class TList>
    const unsigned int
    TypelistBraket<BFS,TList>::num_part() const
  {
    return Length<TList>::value;
  }

  template <class BFS, class TList>
    const unsigned int
    TypelistBraket<BFS,TList>::num_members(unsigned int p) const
    {
    }
  
  /*template <class BFS, class TList>
    bool
    TypelistBraket<BFS, TList>::equiv(const TypelistBraket<BFS,TList>& a) const
  {
    return bfs_ == a.bfs_;
  }*/
  

};

#endif
