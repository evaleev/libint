
#include <string>
#include <hashable.h>
#include <global_macros.h>

#ifndef _libint2_src_bin_libint_oper_h_
#define _libint2_src_bin_libint_oper_h_

namespace libint2 {
  
  /** Permutational symmetries: antisymmetric(anti), symmetric(symm), nonsymmetric (nonsymm),
      some more complicated symmetry (nonstd) */
  typedef struct {
    typedef enum {anti=-1, symm=1, nonsymm=0, nonstd=-2} type;
  } PermutationalSymmetry;
  
  /** OperatorProperties describes various properties of an operator or operator set
      np -- number of particles
      multi -- true if multiplicative
    */
  template <unsigned int NP, bool multi, PermutationalSymmetry::type psymmetry>
    class OperatorProperties {
    public:
      static const unsigned int np = NP;
      static const bool multiplicative = multi;
      static const PermutationalSymmetry::type psymm = psymmetry;
    };

  /** OperSet is the base class for all (sets of) operators.
     OperSet's must be constructable using
     SafePtr<OperSet> or SafePtr<ConstructablePolymorphically>.
  */
  class OperSet : public ConstructablePolymorphically {
    public:
      virtual ~OperSet() {};

      /// Returns full description of the operator
      virtual const std::string& descr() const =0;
      /// Returns short label for the operator
      virtual const std::string& label() const =0;
      /** Returns 1, 0, or -1, if each operator in the set is symmetric, nonsymmetric,
        or antisymmetric with respect to permutation of particles i and j */
      virtual int psymm(int i, int j) const =0;

      /// Number of operators in the set
      virtual const unsigned int num_oper() const =0;
    };
  
  /**
  */
  template <class Props>
    class Oper : public OperSet {
    public:
      typedef Props Properties;
      virtual ~Oper();

      /// Implementation of OperSet::descr()
      const std::string& descr() const;
      /// Implementation of OperSet::label()
      const std::string& label() const;
      /// Implementation of OperSet::psymm()
      int psymm(int i, int j) const;
      
      bool operator==(const Oper&) const;

    protected:
      /// The only declared constructor is only useable by derived classes
      Oper(const std::string& descr, const std::string& label);

    private:
      /// Described name
      const std::string descr_;
      /// short (<20 chars) ID label
      const std::string label_;
      
      /// Must be overloaded by derived classes if Props::psymm == PermutationalSymmetry::nonstd
      virtual int nonstd_psymm(int i, int j) const;
  };

  template <class Props>
    Oper<Props>::Oper(const std::string& descr, const std::string& label) :
    OperSet(), descr_(descr), label_(label)
    {
    }
  
  template <class Props>
    Oper<Props>::~Oper()
    {
    }
  
  template <class Props>
    const std::string&
    Oper<Props>::descr() const
    {
      return descr_;
    }
  
  template <class Props>
    const std::string&
    Oper<Props>::label() const
    {
      return label_;
    }
  
  template <class Props>
    int
    Oper<Props>::psymm(int i, int j) const
    {
      if (i<0 || i>=Props::np)
        throw std::runtime_error("Oper<Props>::psymm(i,j) -- index i out of bounds");
      if (j<0 || j>=Props::np)
        throw std::runtime_error("Oper<Props>::psymm(i,j) -- index j out of bounds");
      if (i == j)
        return 1;
      
      switch(Props::psymm) {
        case PermutationalSymmetry::anti:
        return -1;
        case PermutationalSymmetry::symm:
        return 1;
        case PermutationalSymmetry::nonsymm:
        return 0;
        case PermutationalSymmetry::nonstd:
        return nonstd_psymm(i,j);
      }
    }
  
  template <class Props>
    int
    Oper<Props>::nonstd_psymm(int i, int j) const
    {
      throw ProgrammingError("Props::psymm == nonstd but nonstd_psymm is not overloaded");
    }
    
  template <class Props>
  bool
    Oper<Props>::operator==(const Oper& a) const
    {
      return true;
    }

//////////////////////////////
  
  /** GenOper is a generic operator
  */
  template <class Props>
    class GenSymmOper : public Oper<Props>, public Hashable<unsigned,ComputeKey> {
      public:
      typedef Oper<Props> parent_type;
      /// GenOper is not a set
      typedef GenSymmOper iter_type;
      const unsigned int num_oper() const { return 1; };
      /// Implementation of Hashable::key()
      unsigned key() const { return 0; }
      /// Range of key is [0,1)
      static const unsigned max_key = 1;
      
      GenSymmOper() : Oper<Props>(std::string("General Symmetric Operator"),std::string("GenSymmOper")) {}
      GenSymmOper(const SafePtr<GenSymmOper>&) : Oper<Props>(std::string("General Symmetric Operator"),std::string("GenSymmOper")) {}
      ~GenSymmOper() {}
      
    };

//////////////////////////////
  
  typedef OperatorProperties<2,true,PermutationalSymmetry::symm> MultiplicativeSymm2Body_Props;
  typedef OperatorProperties<2,true,PermutationalSymmetry::nonsymm> MultiplicativeNonsymm2Body_Props;
  
  /** TwoERep is the two-body repulsion operator.
  */
  class TwoERep : public Oper<MultiplicativeSymm2Body_Props>,
                  public Hashable<unsigned,ComputeKey> {
  public:
    typedef Oper<MultiplicativeSymm2Body_Props> parent_type;
    /// TwoERep is not a set
    typedef TwoERep iter_type;
    const unsigned int num_oper() const { return 1; };
    /// Implementation of Hashable::key()
    unsigned key() const { return 0; }
    /// key is in range [0,1)
    static const unsigned max_key = 1;
  
    TwoERep();
    TwoERep(const SafePtr<TwoERep>&);
    TwoERep(const SafePtr<OperSet>&);
    TwoERep(const SafePtr<ConstructablePolymorphically>&);
    TwoERep(const ConstructablePolymorphically&);
    ~TwoERep();
  };
  
//////////////
  
  /** R12_k_G12 is a two-body operator of form r_{12}^k * exp(-\gamma * r_{12}),
      where k is an integer and \gamma is a positive real number.
  */
  template <int K>
  class R12_k_G12 : public Oper<MultiplicativeSymm2Body_Props>,
                    public Hashable<unsigned,ComputeKey> {
  public:
    typedef Oper<MultiplicativeSymm2Body_Props> parent_type;
    /// R12_k_G12 is not a set
    typedef R12_k_G12 iter_type;
    static const int k = K;
    const unsigned int num_oper() const { return 1; }
    /// Implementation of Hashable::key()
    unsigned key() const { return 0; }
    /// key is in range [0,1)
    static const unsigned max_key = 1;
    
    R12_k_G12();
    R12_k_G12(const SafePtr<R12_k_G12>&);
    R12_k_G12(const SafePtr<OperSet>&);
    R12_k_G12(const SafePtr<ConstructablePolymorphically>&);
    R12_k_G12(const ConstructablePolymorphically&);
    ~R12_k_G12();
  };

  namespace R12kG12 {
    std::string label(int K);
    std::string symbol(int K);
  };

  template <int K>
  R12_k_G12<K>::R12_k_G12() :
  parent_type(R12kG12::label(K),R12kG12::symbol(K))
  {
  }
  
  template <int K>
  R12_k_G12<K>::R12_k_G12(const SafePtr<R12_k_G12>& source) :
  parent_type(R12kG12::label(K),R12kG12::symbol(K))
  {
  }
  
  template <int K>
  R12_k_G12<K>::R12_k_G12(const SafePtr<OperSet>& oset) :
  parent_type(R12kG12::label(K),R12kG12::symbol(K))
  {
    const SafePtr<R12_k_G12> oset_cast = dynamic_pointer_cast<R12_k_G12,OperSet>(oset);
    if (oset_cast == 0)
      throw std::runtime_error("R12_k_G12<K>::R12_k_G12(const SafePtr<OperSet>& oset) -- oset is a pointer to an incompatible type");
  }
  
  template <int K>
  R12_k_G12<K>::R12_k_G12(const SafePtr<ConstructablePolymorphically>& oset) :
  parent_type(R12kG12::label(K),R12kG12::symbol(K))
  {
    const SafePtr<R12_k_G12> oset_cast = dynamic_pointer_cast<R12_k_G12,ConstructablePolymorphically>(oset);
    if (oset_cast == 0)
      throw std::runtime_error("R12_k_G12<K>::R12_k_G12(const SafePtr<ConstructablePolymorphically>& oset) -- oset is a pointer to an incompatible type");
  }
  
  template <int K>
  R12_k_G12<K>::R12_k_G12(const ConstructablePolymorphically& oset) :
  parent_type(R12kG12::label(K),R12kG12::symbol(K))
  {
    const R12_k_G12& oset_cast = dynamic_cast<const R12_k_G12&>(oset);
  }

  template <int K>
  R12_k_G12<K>::~R12_k_G12()
  {
  }

//////////////
  
  typedef OperatorProperties<2,false,PermutationalSymmetry::nonsymm> NonmultiplicativeNonsymm2Body_Props;

  /** Ti_G12 is a two-body operator of form [T_i, G12],
      where i is particle index (0 or 1) and G12 is a Gaussian Geminal.
  */
  template <int I>
  class Ti_G12 : public Oper<NonmultiplicativeNonsymm2Body_Props>,
                 public Hashable<unsigned,ComputeKey> {
  public:
    typedef Oper<NonmultiplicativeNonsymm2Body_Props> parent_type;
    /// Ti_G12 is not a set
    typedef Ti_G12 iter_type;
    static const int i = I;
    const unsigned int num_oper() const { return 1; }
    /// Implementation of Hashable::key()
    unsigned key() const { return 0; }
    /// key is in range [0,1)
    static const unsigned max_key = 1;
    
    Ti_G12();
    Ti_G12(const SafePtr<Ti_G12>&);
    Ti_G12(const SafePtr<OperSet>&);
    Ti_G12(const SafePtr<ConstructablePolymorphically>&);
    Ti_G12(const ConstructablePolymorphically&);
    ~Ti_G12();

  };

  namespace TiG12 {
    std::string label(int K);
    std::string symbol(int K);
  };

  template <int I>
  Ti_G12<I>::Ti_G12() :
  parent_type(TiG12::label(I),TiG12::symbol(I))
  {
  }
  
  template <int I>
  Ti_G12<I>::Ti_G12(const SafePtr<Ti_G12>& source) :
  parent_type(TiG12::label(I),TiG12::symbol(I))
  {
  }
  
  template <int I>
  Ti_G12<I>::Ti_G12(const SafePtr<OperSet>& oset) :
  parent_type(TiG12::label(I),TiG12::symbol(I))
  {
    const SafePtr<Ti_G12> oset_cast = dynamic_pointer_cast<Ti_G12,OperSet>(oset);
    if (oset_cast == 0)
      throw std::runtime_error("Ti_G12<I>::Ti_G12(const SafePtr<OperSet>& oset) -- oset is a pointer to an incompatible type");
  }
  
  template <int I>
  Ti_G12<I>::Ti_G12(const SafePtr<ConstructablePolymorphically>& oset) :
  parent_type(TiG12::label(I),TiG12::symbol(I))
  {
    const SafePtr<Ti_G12> oset_cast = dynamic_pointer_cast<Ti_G12,ConstructablePolymorphically>(oset);
    if (oset_cast == 0)
      throw std::runtime_error("Ti_G12<I>::Ti_G12(const SafePtr<ConstructablePolymorphically>& oset) -- oset is a pointer to an incompatible type");
  }
  
  template <int I>
  Ti_G12<I>::Ti_G12(const ConstructablePolymorphically& oset) :
  parent_type(TiG12::label(I),TiG12::symbol(I))
  {
    const Ti_G12& oset_cast = dynamic_cast<const Ti_G12&>(oset);
  }

  template <int I>
  Ti_G12<I>::~Ti_G12()
  {
  }

//////////////

  /** r_1.r_1 x g12 is a result of differentiation of exp( - a r_1^2 - a r_2^2 - c r_{12}^2) geminal .
  */
  class R1dotR1_G12 : public Oper<MultiplicativeNonsymm2Body_Props>,
                      public Hashable<unsigned,ComputeKey> {
  public:
    typedef Oper<MultiplicativeNonsymm2Body_Props> parent_type;
    /// R1dotR1_G12 is not a set
    typedef R1dotR1_G12 iter_type;
    const unsigned int num_oper() const { return 1; };
    /// Implementation of Hashable::key()
    unsigned key() const { return 0; }
    /// key is in range [0,1)
    static const unsigned max_key = 1;
  
    R1dotR1_G12();
    R1dotR1_G12(const SafePtr<R1dotR1_G12>&);
    R1dotR1_G12(const SafePtr<OperSet>&);
    R1dotR1_G12(const SafePtr<ConstructablePolymorphically>&);
    R1dotR1_G12(const ConstructablePolymorphically&);
    ~R1dotR1_G12();
  };

  /** r_1.r_2 x g12 is a result of differentiation of exp( - a r_1^2 - a r_2^2 - c r_{12}^2) geminal .
  */
  class R1dotR2_G12 : public Oper<MultiplicativeSymm2Body_Props>,
                      public Hashable<unsigned,ComputeKey> {
  public:
    typedef Oper<MultiplicativeSymm2Body_Props> parent_type;
    /// R1dotR2_G12 is not a set
    typedef R1dotR2_G12 iter_type;
    const unsigned int num_oper() const { return 1; };
    /// Implementation of Hashable::key()
    unsigned key() const { return 0; }
    /// key is in range [0,1)
    static const unsigned max_key = 1;
  
    R1dotR2_G12();
    R1dotR2_G12(const SafePtr<R1dotR2_G12>&);
    R1dotR2_G12(const SafePtr<OperSet>&);
    R1dotR2_G12(const SafePtr<ConstructablePolymorphically>&);
    R1dotR2_G12(const ConstructablePolymorphically&);
    ~R1dotR2_G12();
  };

};

#endif

