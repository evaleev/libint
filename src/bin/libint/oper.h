
#include <string>

#ifndef _libint2_src_bin_libint_oper_h_
#define _libint2_src_bin_libint_oper_h_

namespace libint2 {
  
  /** OperatorProperties describes various properties of an operator or operator set
      np -- number of particles
      multi -- true if multiplicative
    */
  template <unsigned int NP, bool multi>
    class OperatorProperties {
    public:
      static const unsigned int np = NP;
      static const bool multiplicative = multi;
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

      /// Returns full description of the operator
      const std::string& descr() const;
      /// Returns short label for the operator
      const std::string& label() const;
      
      bool operator==(const Oper&) const;

    protected:
      /// The only declared constructor is only useable by derived classes
      Oper(const std::string& descr, const std::string& label);

    private:
      /// Described name
      const std::string descr_;
      /// short (<20 chars) ID label
      const std::string label_;
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
  bool
    Oper<Props>::operator==(const Oper& a) const
    {
      return true;
    }


//////////////////////////////
  
  typedef OperatorProperties<2,true> Multiplicative2Body_Props;
  
  /** TwoERep is the two-body repulsion operator.
  */
  class TwoERep : public Oper<Multiplicative2Body_Props> {
  public:
    typedef Oper<Multiplicative2Body_Props> parent_type;
    /// TwoERep is not a set
    typedef TwoERep iter_type;
    const unsigned int num_oper() const { return 1; };
  
    TwoERep();
    TwoERep(const SafePtr<TwoERep>&);
    TwoERep(const SafePtr<OperSet>&);
    TwoERep(const SafePtr<ConstructablePolymorphically>&);
    ~TwoERep();

    /** Returns 1, 0, or -1, if the operator is symmetric, nonsymmetric,
        or antisymmetric with respect to permutation of particles i and j */
    int psymm(int i, int j) const;

  private:
    // symmetry W.R.T. permutation of each pair of particles
    // 1 -- symmetric, -1 -- antisymmetric, 0 -- nonsymmetric
    // stored as a lower triangle (diagonal not included)
    static const char psymm_[Properties::np*(Properties::np-1)/2];

  };
  
//////////////
  
  /** R12_k_G12 is a two-body operator of form r_{12}^k * exp(-\gamma * r_{12}),
      where k is an integer and \gamma is a positive real number.
  */
  template <int K>
  class R12_k_G12 : public Oper<Multiplicative2Body_Props> {
  public:
    typedef Oper<Multiplicative2Body_Props> parent_type;
    /// R12_k_G12 is not a set
    typedef R12_k_G12 iter_type;
    static const int k = K;
    const unsigned int num_oper() const { return 1; };
    
    R12_k_G12();
    R12_k_G12(const SafePtr<R12_k_G12>&);
    R12_k_G12(const SafePtr<OperSet>&);
    R12_k_G12(const SafePtr<ConstructablePolymorphically>&);
    ~R12_k_G12();

    /** Returns 1, 0, or -1, if the operator is symmetric, nonsymmetric,
        or antisymmetric with respect to permutation of particles i and j */
    int psymm(int i, int j) const;

  private:
    // symmetry W.R.T. permutation of each pair of particles
    // 1 -- symmetric, -1 -- antisymmetric, 0 -- nonsymmetric
    // stored as a lower triangle (diagonal not included)
    static const int npair = Properties::np*(Properties::np-1)/2;
    static const vector<char> psymm_;

  };
  
  template <int K> const vector<char> R12_k_G12<K>::psymm_(npair,1);

  template <int K>
  R12_k_G12<K>::R12_k_G12() :
  parent_type("R12^k * G12","R12_k_G12")
  {
  }
  
  template <int K>
  R12_k_G12<K>::R12_k_G12(const SafePtr<R12_k_G12>& source) :
  parent_type("R12^k * G12","R12_k_G12")
  {
  }
  
  template <int K>
  R12_k_G12<K>::R12_k_G12(const SafePtr<OperSet>& oset) :
  parent_type("R12^k * G12","R12_k_G12")
  {
    const SafePtr<R12_k_G12> oset_cast = dynamic_pointer_cast<R12_k_G12,OperSet>(oset);
    if (oset_cast == 0)
      throw std::runtime_error("R12_k_G12<K>::R12_k_G12(const SafePtr<OperSet>& oset) -- oset is a pointer to an incompatible type");
  }
  
  template <int K>
  R12_k_G12<K>::R12_k_G12(const SafePtr<ConstructablePolymorphically>& oset) :
  parent_type("R12^k * G12","R12_k_G12")
  {
    const SafePtr<R12_k_G12> oset_cast = dynamic_pointer_cast<R12_k_G12,ConstructablePolymorphically>(oset);
    if (oset_cast == 0)
      throw std::runtime_error("R12_k_G12<K>::R12_k_G12(const SafePtr<ConstructablePolymorphically>& oset) -- oset is a pointer to an incompatible type");
  }
  
  template <int K>
  R12_k_G12<K>::~R12_k_G12()
  {
  }
  
  template <int K>
  int
  R12_k_G12<K>::psymm(int i, int j) const
  {
    if (i<0 || i>=Properties::np)
      throw std::runtime_error("R12_k_G12<K>::psymm(i,j) -- index i out of bounds");
    if (j<0 || j>=Properties::np)
      throw std::runtime_error("R12_k_G12<K>::psymm(i,j) -- index j out of bounds");
    if (i == j)
      return 1;
    int ii = (i > j) ? i : j;
    int jj = (i > j) ? j : i;
    return psymm_[ii*(ii-1)/2 + jj];
  }

//////////////
  
  typedef OperatorProperties<2,false> Nonmultiplicative2Body_Props;

  /** Ti_G12 is a two-body operator of form [T_i, G12],
      where i is particle index (0 or 1) and G12 is a Gaussian Geminal.
  */
  template <int I>
  class Ti_G12 : public Oper<Nonmultiplicative2Body_Props> {
  public:
    typedef Oper<Nonmultiplicative2Body_Props> parent_type;
    /// Ti_G12 is not a set
    typedef Ti_G12 iter_type;
    static const int i = I;
    const unsigned int num_oper() const { return 1; };
    
    Ti_G12();
    Ti_G12(const SafePtr<Ti_G12>&);
    Ti_G12(const SafePtr<OperSet>&);
    Ti_G12(const SafePtr<ConstructablePolymorphically>&);
    ~Ti_G12();

    /** Returns 1, 0, or -1, if the operator is symmetric, nonsymmetric,
        or antisymmetric with respect to permutation of particles i and j */
    int psymm(int i, int j) const;

  private:
    // symmetry W.R.T. permutation of each pair of particles
    // 1 -- symmetric, -1 -- antisymmetric, 0 -- nonsymmetric
    // stored as a lower triangle (diagonal not included)
    static const int npair = Properties::np*(Properties::np-1)/2;
    static const vector<char> psymm_;
  };
  
  template <int I> const vector<char> Ti_G12<I>::psymm_(npair,0);

  template <int I>
  Ti_G12<I>::Ti_G12() :
  parent_type("[T_i,G12]","Ti_G12")
  {
  }
  
  template <int I>
  Ti_G12<I>::Ti_G12(const SafePtr<Ti_G12>& source) :
  parent_type("[T_i,G12]","Ti_G12")
  {
  }
  
  template <int I>
  Ti_G12<I>::Ti_G12(const SafePtr<OperSet>& oset) :
  parent_type("[T_i,G12]","Ti_G12")
  {
    const SafePtr<Ti_G12> oset_cast = dynamic_pointer_cast<Ti_G12,OperSet>(oset);
    if (oset_cast == 0)
      throw std::runtime_error("Ti_G12<I>::Ti_G12(const SafePtr<OperSet>& oset) -- oset is a pointer to an incompatible type");
  }
  
  template <int I>
  Ti_G12<I>::Ti_G12(const SafePtr<ConstructablePolymorphically>& oset) :
  parent_type("[T_i,G12]","Ti_G12")
  {
    const SafePtr<Ti_G12> oset_cast = dynamic_pointer_cast<Ti_G12,ConstructablePolymorphically>(oset);
    if (oset_cast == 0)
      throw std::runtime_error("Ti_G12<I>::Ti_G12(const SafePtr<ConstructablePolymorphically>& oset) -- oset is a pointer to an incompatible type");
  }
  
  template <int I>
  Ti_G12<I>::~Ti_G12()
  {
  }
  
  template <int I>
  int
  Ti_G12<I>::psymm(int i, int j) const
  {
    if (i<0 || i>=Properties::np)
      throw std::runtime_error("Ti_G12<I>::psymm(i,j) -- index i out of bounds");
    if (j<0 || j>=Properties::np)
      throw std::runtime_error("Ti_G12<I>::psymm(i,j) -- index j out of bounds");
    if (i == j)
      return 1;
    int ii = (i > j) ? i : j;
    int jj = (i > j) ? j : i;
    return psymm_[ii*(ii-1)/2 + jj];
  }


};

#endif

