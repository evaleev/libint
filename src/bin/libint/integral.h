
#include <smart_ptr.h>
#include <rr.h>
#include <typelist.h>
#include <iter.h>
#include <quanta.h>
#include <equiv.h>

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

    /// Obtain pointers to ith BasisFunctionSet for particle p in bra
    virtual const SafePtr<BasisFunctionSet> bra(unsigned int p, unsigned int i) const =0;
    /// Obtain pointers to ith BasisFunctionSet for particle p in ket
    virtual const SafePtr<BasisFunctionSet> ket(unsigned int p, unsigned int i) const =0;

  };

  /**
     GenIntegralSet is a set of integrals over functions derived from BFS.

     Oper is an operator or a set of operators. Oper must be derived from OperSet.
     BraSetType and KetSetType are types describing sets of functions.
     An example of a class that can be used as BraSetType and KetSetType
     is VectorBraket.
     AuxQuanta describes auxiliary quantum numbers. AuxQuanta should be derived from QuantumSet.
  */
  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta = NullQuantumSet> class GenIntegralSet :
    public IntegralSet<BFS>, public DGVertex
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
    
      /** No constructors are public since this is a singleton-like quantity.
          Instead, access is provided through Instance().
          
          Derived classes do not have to be Singletons -- hence protected constructors are provided.
      */
      virtual ~GenIntegralSet();

      /// Returns a pointer to a unique instance, a la Singleton
      static const SafePtr<GenIntegralSet> Instance(const Oper& oper, const BraSetType& bra, const KetSetType& ket, const AuxQuanta& aux);
      
      /// Comparison operator
      virtual bool operator==(const GenIntegralSet&) const;
      /// Specialization of DGVertex's equiv
      bool equiv(const SafePtr<DGVertex>& v) const { return PtrComp::equiv(this,v); }
      /// Specialization of DGVertex's print
      virtual void print(std::ostream& os = std::cout) const;

      /// Obtain BFsets members
      const SafePtr<BFS> bra(unsigned int p, unsigned int i) const;
      const SafePtr<BFS> ket(unsigned int p, unsigned int i) const;
      
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

    protected:
      // Basic Integral constructor. It is protected so that derived classes don't have to behave like singletons
      GenIntegralSet(const Oper& oper, const BraSetType& bra, const KetSetType& ket, const AuxQuanta& aux);
      
      BraSetType bra_;
      KetSetType ket_;

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

      // Unique instances of GenIntegralSet are placed here and obtained through Instance()
      static vector < SafePtr<GenIntegralSet> > stack_;
      
      // The operator needs to be a real object rather than real type to be able to construct a SubIterator, etc.
      SafePtr<Oper> O_;
      // Same for AuxQuanta
      SafePtr<AuxQuanta> aux_;

    };

  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    vector < SafePtr< GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta> > > GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta>::stack_(0);
  
  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta>::GenIntegralSet(const Oper& oper, const BraSetType& bra, const KetSetType& ket, const AuxQuanta& aux) :
    O_(SafePtr<Oper>(new Oper(oper))), bra_(bra), ket_(ket), aux_(SafePtr<AuxQuanta>(new AuxQuanta(aux)))
    {
      if (Oper::Properties::np != bra.num_part())
        throw std::runtime_error("GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta>::GenIntegralSet(bra,ket) -- number of particles in bra doesn't match that in the operator");
      if (Oper::Properties::np != ket.num_part())
        throw std::runtime_error("GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta>::GenIntegralSet(bra,ket) -- number of particles in ket doesn't match that in the operator");
    }

  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta>::~GenIntegralSet()
    {
    }

  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const SafePtr< GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta> >
    GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta>::Instance(const Oper& oper, const BraSetType& bra, const KetSetType& ket, const AuxQuanta& aux)
    {
      SafePtr<this_type> this_int(new this_type(oper,bra,ket,aux));
      int stack_size = stack_.size();
      for(int i=0; i<stack_size; i++) {
        if (PtrComp::equiv(this_int,stack_[i])) {
          this_int.reset();
          return stack_[i];
        }
      }
      stack_.push_back(this_int);
      return this_int;
    }
  
  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const SafePtr<BFS>
    GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta>::bra(unsigned int p, unsigned int i) const
    {
      return bra_.member(p,i);
    }
    
  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const SafePtr<BFS>
    GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta>::ket(unsigned int p, unsigned int i) const
    {
      return ket_.member(p,i);
    }

  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const typename GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta>::BraType&
    GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta>::bra() const
    {
      return bra_;
    }

  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const typename GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta>::KetType&
    GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta>::ket() const
    {
      return ket_;
    }

  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const SafePtr<Oper>
    GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta>::oper() const
    {
      return O_;
    }

  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const SafePtr<AuxQuanta>
    GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta>::aux() const
    {
      return aux_;
    }

  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    bool
    GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta>::operator==(const this_type& a) const
    {
      bool oper_equiv = PtrEquiv<Oper>::equiv(O_,a.O_);
      bool bra_equiv = PtrEquiv<BraSetType>::equiv(bra_,a.bra_);
      bool ket_equiv = PtrEquiv<KetSetType>::equiv(ket_,a.ket_);
      bool aux_equiv = PtrEquiv<AuxQuanta>::equiv(aux_,a.aux_);
      return oper_equiv && bra_equiv && ket_equiv && aux_equiv;
    }

  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    void
    GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta>::print(std::ostream& os) const
    {
      os << "GenIntegralSet: < ";
      for(int p=0; p<Oper::Properties::np; p++) {
        unsigned int nbra = bra_.num_members(p);
        for(unsigned int i=0; i<nbra; i++)
          os << bra_.member(p,i)->label() << "(" << p << ") ";
      }
      os << "| " << O_->label() << " | ";
      for(int p=0; p<Oper::Properties::np; p++) {
        unsigned int nket = ket_.num_members(p);
        for(unsigned int i=0; i<nket; i++)
          os << ket_.member(p,i)->label() << "(" << p << ") ";
        os << "> ^ { " << aux_->label() << " }" << endl;
      }
    }
  
  /** VectorBraket is a std::vector-based type that can be used as a BraSetType or a KetSetType parameter
      to construct an instance of GenIntegralSet
  */
  template <class BFS> class VectorBraket {

  public:
    typedef vector< SafePtr<BFS> > BFSVector;
    typedef vector< BFSVector > BFSMatrix;
    typedef VectorBraket<typename BFS::iter_type> iter_type;
    typedef struct{} parent_type;

    /** This one is a very dangerous constructor -- do not to use it if at all possible.
      Provided only for compatibility for generic subiterator algorithms */
    VectorBraket();
    VectorBraket(const BFSMatrix&);
    VectorBraket(const VectorBraket&);
    ~VectorBraket() throw();

    /// Comparison function
    bool operator==(const VectorBraket&) const;
    /// Returns base pointer to the i-th function for particle p
    //const SafePtr<ConstructablePolymorphically> member(unsigned int p, unsigned int i) const;
    /// Returns pointer to the i-th function for particle p
    const SafePtr<BFS> member(unsigned int p, unsigned int i) const;
    /// Returns pointer to the SubIterator for i-th BFS of particle p
    SubIterator* member_subiter(unsigned int p, unsigned int i) const;
    /// Sets i-th function for particle p
    void set_member(const BFS&, unsigned int p, unsigned int i);
    /// Sets i-th function for particle p
    void set_member(const SafePtr<BFS>&, unsigned int p, unsigned int i);
    /// Sets i-th function for particle p (does a dynamic cast inside)
    void set_member(const SafePtr<ConstructablePolymorphically>&, unsigned int p, unsigned int i);
    /// Returns the number of BFS for particle p
    const unsigned int num_members(unsigned int p) const;
    /// Returns the number of particles
    const unsigned int num_part() const;

  private:
    
    BFSMatrix bfs_;

  };

  template <class BFS>
    VectorBraket<BFS>::VectorBraket() :
    bfs_(0)
    {
    }
  
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

  /*
  template <class BFS>
    const SafePtr<ConstructablePolymorphically>
    VectorBraket<BFS>::member(unsigned int p, unsigned int i) const
    {
      const SafePtr<ConstructablePolymorphically> ptr = dynamic_pointer_cast<ConstructablePolymorphically,BFS>(bfs_.at(p).at(i));
      return ptr;
    }
  */

  template <class BFS>
    const SafePtr<BFS>
    VectorBraket<BFS>::member(unsigned int p, unsigned int i) const
    {
      return bfs_.at(p).at(i);
    }

  template <class BFS>
    SubIterator*
    VectorBraket<BFS>::member_subiter(unsigned int p, unsigned int i) const
    {
      return static_cast<SubIterator*>(new SubIteratorBase<BFS>( member(p,i) ) );
    }
  
  template <class BFS>
    void
    VectorBraket<BFS>::set_member(const BFS& bfs, unsigned int p, unsigned int i)
    {
      if (p >= bfs_.size())
        bfs_.resize(p+1);
      if (i >= bfs_[p].size())
        bfs_[p].resize(i+1);
      SafePtr<BFS> bfs_ptr(new BFS(bfs));
      bfs_[p][i] = bfs_ptr;
    }

  template <class BFS>
    void
    VectorBraket<BFS>::set_member(const SafePtr<ConstructablePolymorphically>& bfs, unsigned int p, unsigned int i)
    {
      // WARNING : can be VERY dangerous
      // try constructing BFS from bfs.
      SafePtr<BFS> bfs_cast(new BFS(bfs));
      
      if (p >= bfs_.size())
        bfs_.resize(p+1);
      if (i >= bfs_[p].size())
        bfs_[p].resize(i+1);
      bfs_[p][i] = bfs_cast;
    }
  
  template <class BFS>
    void
    VectorBraket<BFS>::set_member(const SafePtr<BFS>& bfs, unsigned int p, unsigned int i)
    {
      if (p >= bfs_.size())
        bfs_.resize(p+1);
      if (i >= bfs_[p].size())
        bfs_[p].resize(i+1);
      bfs_[p][i] = bfs;
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
    VectorBraket<BFS>::operator==(const VectorBraket<BFS>& a) const
    {
      if (bfs_.size() != a.bfs_.size())
        return false;

      // compare each row
      const int size1 = bfs_.size();
      for(int i=0; i<size1; i++) {
        if (bfs_[i].size() != a.bfs_[i].size())
          return false;

        // compare each element
        const int size2 = bfs_[i].size();
        for(int j=0; j<size2; j++)
          if (!PtrEquiv<BFS>::equiv(bfs_[i][j],a.bfs_[i][j]))
            return false;
      }
      return true;
    }

  /**
     Most basic type -- TwoPRep_11_11 --
     has one bfs for each particle in bra and ket.
     Note that GenIntegralSet is initialized with an abstract type libint2::BFSet,
     from which BFS derives.
  */
  template <class BFS> class TwoPRep_11_11 :
    public GenIntegralSet<TwoERep, IncableBFSet, VectorBraket<BFS>, VectorBraket<BFS>, QuantumNumbers<unsigned int,1> >
    {
    public:
      typedef VectorBraket<BFS> BraType;
      typedef VectorBraket<BFS> KetType;
      typedef QuantumNumbers<unsigned int,1> AuxIndexType;
      typedef TwoPRep_11_11 this_type;
      
      /// TwoPRep_11_11 is a set of these subobjects
      typedef TwoPRep_11_11<typename BFS::iter_type> iter_type;
      /// This is the immediate parent
      typedef GenIntegralSet<TwoERep, IncableBFSet, VectorBraket<BFS>, VectorBraket<BFS>, AuxIndexType > parent_type;
      /// This class provides comparison operations on pointers
      typedef PtrEquiv<this_type> PtrComp;

      /* This "constructor" takes basis function sets, in Mulliken ordering.
         Returns a pointer to a unique instance, a la Singleton
      */
      static const SafePtr<TwoPRep_11_11> Instance(const BFS& bra0, const BFS& ket0, const BFS& bra1, const BFS& ket1, unsigned int m);
      /// Returns a pointer to a unique instance, a la Singleton
      static const SafePtr<TwoPRep_11_11> Instance(const VectorBraket<BFS>& bra, const VectorBraket<BFS>& ket, const AuxIndexType& aux);
      
      unsigned int m() const { return parent_type::aux()->elem(0); };
      
      /// Comparison operator
      bool operator==(const this_type&) const;
      /// Specialization of GenIntegralSet's and DGVertex's print
      void print(std::ostream& os = std::cout) const;

    private:
      // This constructor is also private and not implemented since all Integral's are Singletons. Use Instance instead.
      TwoPRep_11_11(const VectorBraket<BFS>& bra, const VectorBraket<BFS>& ket, const AuxIndexType& aux);

      /// stack_ of pointers to objects used to check whether an object already exists
      static vector< SafePtr<TwoPRep_11_11> > stack_;

    };

  template <class BFS>
    vector< SafePtr< TwoPRep_11_11<BFS> > > TwoPRep_11_11<BFS>::stack_(0);

  template <class BFS>
    TwoPRep_11_11<BFS>::TwoPRep_11_11(const VectorBraket<BFS>& bra, const VectorBraket<BFS>& ket,  const AuxIndexType& aux) :
    GenIntegralSet<TwoERep, IncableBFSet, BraType, KetType, AuxIndexType>(TwoERep(),bra, ket, aux)
    {
      if (bra.num_members(0) != 1)
        throw std::runtime_error("TwoPRep_11_11<BFS>::TwoPRep_11_11(bra,ket) -- number of BFSs in bra for particle 0 must be 1");
      if (bra.num_members(1) != 1)
        throw std::runtime_error("TwoPRep_11_11<BFS>::TwoPRep_11_11(bra,ket) -- number of BFSs in bra for particle 1 must be 1");
      if (ket.num_members(0) != 1)
        throw std::runtime_error("TwoPRep_11_11<BFS>::TwoPRep_11_11(bra,ket) -- number of BFSs in ket for particle 0 must be 1");
      if (ket.num_members(1) != 1)
        throw std::runtime_error("TwoPRep_11_11<BFS>::TwoPRep_11_11(bra,ket) -- number of BFSs in ket for particle 1 must be 1");
    }

  template <class BFS>
    const SafePtr< TwoPRep_11_11<BFS> >
    TwoPRep_11_11<BFS>::Instance(const VectorBraket<BFS>& bra, const VectorBraket<BFS>& ket, const AuxIndexType& aux)
    {
      SafePtr<TwoPRep_11_11> this_int(new TwoPRep_11_11<BFS>(bra,ket,aux));
      int stack_size = stack_.size();
      for(int i=0; i<stack_size; i++) {
        if (PtrComp::equiv(this_int,stack_[i])) {
          this_int.reset();
          return stack_[i];
        }
      }
      stack_.push_back(this_int);
      return this_int;
    }

  template <class BFS>
    const SafePtr< TwoPRep_11_11<BFS> >
    TwoPRep_11_11<BFS>::Instance(const BFS& bra0, const BFS& ket0, const BFS& bra1, const BFS& ket1, unsigned int m)
    {
      typedef SafePtr<BFS> BFSPtr;
      BFSPtr bra0_ptr(new BFS(bra0));
      BFSPtr bra1_ptr(new BFS(bra1));
      BFSPtr ket0_ptr(new BFS(ket0));
      BFSPtr ket1_ptr(new BFS(ket1));
      vector<BFSPtr> vbra0;  vbra0.push_back(bra0_ptr);
      vector<BFSPtr> vbra1;  vbra1.push_back(bra1_ptr);
      vector<BFSPtr> vket0;  vket0.push_back(ket0_ptr);
      vector<BFSPtr> vket1;  vket1.push_back(ket1_ptr);
      vector< vector<BFSPtr> > vvbra;  vvbra.push_back(vbra0);  vvbra.push_back(vbra1);
      vector< vector<BFSPtr> > vvket;  vvket.push_back(vket0);  vvket.push_back(vket1);
      VectorBraket<BFS> bra(vvbra);
      VectorBraket<BFS> ket(vvket);
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
    TwoPRep_11_11<BFS>::print(std::ostream& os) const
    {
      os << "TwoPRep_11_11: (" << parent_type::bra_.member(0,0)->label() << " " << parent_type::ket_.member(0,0)->label()
         << " | " << parent_type::bra_.member(1,0)->label() << " " << parent_type::ket_.member(1,0)->label() << ")^{" << m() <<"}" << endl;
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
