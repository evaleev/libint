
#include <typelist.h>

#ifndef _libint2_src_bin_libint_integral_h_
#define _libint2_src_bin_libint_integral_h_

namespace libint2 {

  /**
     O is an operator
     TList is a typelist TList1 that consists 2 second-level typelists: TList2 for bra and TList3 for ket
     TList2 and TList3 consist of O::np types which specify basis functions for each particle in bra/ket
  */
  
  // Tags for first and second typelists
  typedef Int2Type<1> Tag1;
  typedef Int2Type<2> Tag2;
  typedef Int2Type<2> Tag3;

#if 0
  template < class O > struct IntegralImpl
  {
    BFSet* bra_[O::np];
    BFSet* ket_[O::np];
  };

  // General Integral declaration takes Operator, list of bra function types and list of key functions types
  template < class O, class BraList, class KetList > class Integral;
  
  // unrolling BraList
  template < class O, class T, class U, class KetList>
    class Integral< O, PTYPELIST_2(Tag1,T,U), KetList> :
    public Integral<O, T, KetList>,
    public Integral<O, U, KetList>,
    virtual public IntegralImpl<O>
    {
    public:
      typedef PTYPELIST_2(Tag1,T,U) BraList;

      BFSet* bra[O::np];
      BFSet* ket[O::np];

      Integral(const T& bra_fn)

      ~Integral()
        {
          // This must be a list for bra and ket
          if (Length<BraKetList>::value != 2)
            assert(false);
        }
    };

  // KetList
  template < class O, class T>
    class Integral< O, PTYPELIST_2(Tag1,T,NullType<Tag1>) > :
    public Integral<O, T>
    {
    public:
      typedef PTYPELIST_2(Tag1,T,NullType<Tag1>) KetList;
    };

  // Bra list
  template < class O, class T, class U >
    class Integral< O, PTYPELIST_2(Tag2,T,U) > :
    public Integral<O, T>,
    public Integral<O, U>
    {
    public:
      typedef PTYPELIST_2(Tag2,T,U) BraFuncList;

      T* 
      T bra[O::np];

      ~Integral()
        {
          // This must be a list for bra and ket functions
          if (Length<BraFuncList>::value != O::np)
            assert(false);
        }
    };

  // specialization only exists when TList is tagged with Tag1
  template < class O, class T>
    class Integral< O, PTYPELIST_2(Tag2,T,NullType<Tag2>) > :
    public Integral<O, T>
    {
    public:
      typedef PTYPELIST_2(Tag2,T,NullType<Tag2>) KetFuncList;

      T ket[O::np];

      ~Integral()
        {
          // This must be a list for ket functions
          if (Length<KetFuncList>::value != O::np)
            assert(false);
        }
    };

#endif

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
  template <class BasisFunctionSet> class IntegralSetIter;
  template <class BasisFunctionSet> class IntegralSet {

  public:
    virtual ~IntegralSet() {};

    /// Equivalence operator
    virtual bool equiv(const IntegralSet*) const =0;

    /// Obtain pointers to ith BasisFunctionSet for particle p in bra
    virtual const BasisFunctionSet* bra(unsigned int p, unsigned int i) const =0;
    /// Obtain pointers to ith BasisFunctionSet for particle p in ket
    virtual const BasisFunctionSet* ket(unsigned int p, unsigned int i) const =0;

  };

  /**
     GenIntegralSet uses functions derived from BFS
  */
  template <class Oper, class BFS, class BraSetType, class KetSetType> class GenIntegralSet :
    public IntegralSet<BFS>
    {
    public:
      /** No constructors are public since this is a singleton-like quantity.
          Instead, access is provided through derived class's Instance().
      */
      virtual ~GenIntegralSet();
      
      /// Equivalence operator
      virtual bool equiv(const IntegralSet<BFS>*) const;
        
      /// Obtain BFsets
      const BFS* bra(unsigned int p, unsigned int i) const;
      const BFS* ket(unsigned int p, unsigned int i) const;


    protected:
      // Basic Integral constructor
      GenIntegralSet(const BraSetType& bra, const KetSetType& ket);
      
      BraSetType bra_;
      KetSetType ket_;

    private:
      typedef Oper OperatorType;

      //
      // All integrals are Singletons by nature, therefore they must be treated as such
      // 1) No public constructors are provided
      // 2) protected members are provided to implement Singleton-type functionality
      //
      GenIntegralSet();
      GenIntegralSet(const GenIntegralSet&);
      // Copy is not permitted
      GenIntegralSet& operator=(const GenIntegralSet& source);

    };


  template <class Oper, class BFS, class BraSetType, class KetSetType>
    GenIntegralSet<Oper,BFS,BraSetType,KetSetType>::GenIntegralSet(const BraSetType& bra, const KetSetType& ket) :
    bra_(bra), ket_(ket)
    {
      if (Oper::np != bra.num_part())
        throw std::runtime_error("TwoERep_11_11<BFSet>::TwoERep_11_11(bra,ket) -- number of particles in bra doesn't match that in the operator");
      if (Oper::np != ket.num_part())
        throw std::runtime_error("TwoERep_11_11<BFSet>::TwoERep_11_11(bra,ket) -- number of particles in ket doesn't match that in the operator");
    }

  template <class Oper, class BFS, class BraSetType, class KetSetType>
    GenIntegralSet<Oper,BFS,BraSetType,KetSetType>::~GenIntegralSet()
    {
    }

  template <class Oper, class BFS, class BraSetType, class KetSetType>
    const BFS*
    GenIntegralSet<Oper,BFS,BraSetType,KetSetType>::bra(unsigned int p, unsigned int i) const
    {
      return bra_.member(p,i);
    }
    
  template <class Oper, class BFS, class BraSetType, class KetSetType>
    const BFS*
    GenIntegralSet<Oper,BFS,BraSetType,KetSetType>::ket(unsigned int p, unsigned int i) const
    {
      return ket_.member(p,i);
    }
    
  template <class Oper, class BFS, class BraSetType, class KetSetType>
    bool
    GenIntegralSet<Oper,BFS,BraSetType,KetSetType>::equiv(const IntegralSet<BFS>* a) const
    {
      const GenIntegralSet<Oper,BFS,BraSetType,KetSetType>* a_cast = static_cast< const GenIntegralSet<Oper,BFS,BraSetType,KetSetType>* >(a);
      if (a_cast == 0)
        assert(false);
      bool bra_equiv = bra_.equiv(a_cast->bra_);
      bool ket_equiv = ket_.equiv(a_cast->ket_);
      return bra_equiv && ket_equiv;
    }

  /** VectorBraket is a std::vector-based type that can be used as a BraSetType or a KetSetType parameter
      to construct an instance of GenIntegralSet
  */
  template <class BFS> class VectorBraket {

  public:
    typedef vector<BFS> BFSVector;
    typedef vector< BFSVector > BFSMatrix;

    VectorBraket(const BFSMatrix&);
    VectorBraket(const VectorBraket&);
    ~VectorBraket() throw();

    /// Comparison function
    bool equiv(const VectorBraket&) const;
    /// Returns pointer to the i-th function for particle p
    const BFS* member(unsigned int p, unsigned int i) const;
    /// Returns the number of BFS for particle p
    const unsigned int num_members(unsigned int p) const;
    /// Returns the number of particles
    const unsigned int num_part() const;

  private:
    
    BFSMatrix bfs_;

  };

  template <class BFS>
    VectorBraket<BFS>::VectorBraket(const BFSMatrix& bfs) :
    bfs_(bfs)
    {
    }

  template <class BFS>
    VectorBraket<BFS>::VectorBraket(const VectorBraket& a) :
    bfs_(a.bfs_)
    {
    }

  template <class BFS>
    VectorBraket<BFS>::~VectorBraket() throw()
    {
    }

  template <class BFS>
    const BFS*
    VectorBraket<BFS>::member(unsigned int p, unsigned int i) const
    {
      return &bfs_.at(p).at(i);
    }

  template <class BFS>
    const unsigned int
    VectorBraket<BFS>::num_members(unsigned int p) const
    {
      return bfs_.at(p).size();
    }

  template <class BFS>
    const unsigned int
    VectorBraket<BFS>::num_part() const
    {
      return bfs_.size();
    }

  template <class BFS>
    bool
    VectorBraket<BFS>::equiv(const VectorBraket<BFS>& a) const
    {
      return bfs_ == a.bfs_;
    }

  /**
     Most basis type: TwoERep_11_11
     has one bfs for each particle in bra and ket
  */
  template <class BFS> class TwoERep_11_11 :
    /// Note that GenIntegralSet is initialized with an abstract type libint2::BFSet
    public GenIntegralSet<TwoERep, BFS, VectorBraket<BFS>, VectorBraket<BFS> >,
    public DGVertex
    {
    public:
      typedef VectorBraket<BFS> BraType;
      typedef VectorBraket<BFS> KetType;

      /* This "constructor" takes basis function sets, in Mulliken ordering.
         Returns a pointer to a unique instance, a la Singleton
      */
      static TwoERep_11_11* Instance(const BFS& bra0, const BFS& ket0, const BFS& bra1, const BFS& ket1, unsigned int m);
      
      unsigned int m() const { return m_; };
      
      /// Specialization of DGVertex's equiv
      bool equiv(const DGVertex*) const;

      void print(std::ostream& os = std::cout) const;
      

    private:
      unsigned int m_;  // auxiliary index

      // Default and copy constructors are not allowed
      TwoERep_11_11();
      TwoERep_11_11(const TwoERep_11_11&);      
      // This constructor is also private and not implemented since all Integral's are Singletons. Use Instance instead.
      TwoERep_11_11(const VectorBraket<BFS>& bra, const VectorBraket<BFS>& ket, unsigned int m);

      /// Returns a pointer to a unique instance, a la Singleton
      static TwoERep_11_11* Instance(const VectorBraket<BFS>& bra, const VectorBraket<BFS>& ket, unsigned int m);
      /// stack_ of pointers to objects used to check whether an object already exists
      static vector< TwoERep_11_11* > stack_;

    };

  template <class BFS>
    vector< TwoERep_11_11<BFS>* > TwoERep_11_11<BFS>::stack_(0);

  template <class BFSet>
    TwoERep_11_11<BFSet>::TwoERep_11_11(const VectorBraket<BFSet>& bra, const VectorBraket<BFSet>& ket, unsigned int m) :
    GenIntegralSet<TwoERep, BFSet, BraType, KetType>(bra, ket), DGVertex(), m_(m)
    {
      if (bra.num_members(0) != 1)
        throw std::runtime_error("TwoERep_11_11<BFSet>::TwoERep_11_11(bra,ket) -- number of BFSets in bra for particle 0 must be 1");
      if (bra.num_members(1) != 1)
        throw std::runtime_error("TwoERep_11_11<BFSet>::TwoERep_11_11(bra,ket) -- number of BFSets in bra for particle 1 must be 1");
      if (ket.num_members(0) != 1)
        throw std::runtime_error("TwoERep_11_11<BFSet>::TwoERep_11_11(bra,ket) -- number of BFSets in ket for particle 0 must be 1");
      if (ket.num_members(1) != 1)
        throw std::runtime_error("TwoERep_11_11<BFSet>::TwoERep_11_11(bra,ket) -- number of BFSets in ket for particle 1 must be 1");
    }

  template <class BFSet>
    TwoERep_11_11<BFSet>*
    TwoERep_11_11<BFSet>::Instance(const VectorBraket<BFSet>& bra, const VectorBraket<BFSet>& ket, unsigned int m)
    {
      TwoERep_11_11* const this_int = new TwoERep_11_11<BFSet>(bra,ket,m);
      int stack_size = stack_.size();
      for(int i=0; i<stack_size; i++) {
        if (this_int->equiv(stack_[i])) {
          delete this_int;
          return stack_[i];
        }
      }
      stack_.push_back(this_int);
      return this_int;
    }

  template <class BFSet>
    TwoERep_11_11<BFSet>*
    TwoERep_11_11<BFSet>::Instance(const BFSet& bra0, const BFSet& ket0, const BFSet& bra1, const BFSet& ket1, unsigned int m)
    {
      vector<BFSet> vbra0;  vbra0.push_back(bra0);
      vector<BFSet> vbra1;  vbra1.push_back(bra1);
      vector<BFSet> vket0;  vket0.push_back(ket0);
      vector<BFSet> vket1;  vket1.push_back(ket1);
      vector< vector<BFSet> > vvbra;  vvbra.push_back(vbra0);  vvbra.push_back(vbra1);
      vector< vector<BFSet> > vvket;  vvket.push_back(vket0);  vvket.push_back(vket1);
      VectorBraket<BFSet> bra(vvbra);
      VectorBraket<BFSet> ket(vvket);
      return Instance(bra,ket,m);
    }
  
  template <class BFSet>
    bool
    TwoERep_11_11<BFSet>::equiv(const DGVertex* a) const
    {
      // check the type first
      const TwoERep_11_11<BFSet>* a_cast = dynamic_cast< const TwoERep_11_11<BFSet>* >(a);
      if (!a_cast)
        return false;

      bool result = GenIntegralSet<TwoERep, BFSet, BraType, KetType>::equiv(a_cast) && (m_ == a_cast->m_);
      return result;
    }

  template <class BFSet>
    void
    TwoERep_11_11<BFSet>::print(std::ostream& os) const
    {
      os << "TwoERep_11_11: (" << bra_.member(0,0)->label() << " " << ket_.member(0,0)->label()
         << " | " << bra_.member(1,0)->label() << " " << ket_.member(1,0)->label() << ")^{" << m_ <<"}" << endl;
    };



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
