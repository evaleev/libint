

#ifndef _libint2_src_bin_libint_rr_h_
#define _libint2_src_bin_libint_rr_h_

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

using namespace std;


namespace libint2 {

  const int Libint2_DefaultVectorLength = 64;

  class GaussianShell {
    int l_;
    
    public:
    GaussianShell(int l) {l_ = l;};
    ~GaussianShell() {};
    
  };
  
  class CartesianGaussian : public GaussianShell {

    // cartesian exponents
    int n_[3];

    public:
    CartesianGaussian(int nx, int ny, int nz);
    ~CartesianGaussian() {};

  };

  class Operator {

    // Described name
    const std::string descr_;
    // short (<20 chars) ID label
    const std::string id_;
    
    // number of particles > 0
    char np_;
    // symmetry W.R.T. permutation of each pair of particles
    // 1 -- symmetric, -1 -- antisymmetric, 0 -- nonsymmetric
    // stored as a lower triangle (diagonal not included)
    vector<char> psymm_;
    
    public:
    Operator(const std::string& descr, const std::string& id, char np,
             const vector<char>& psymm);
    ~Operator();

    /// Returns full description of the operator
    const std::string& descr() const;
    /// Returns short label for the operator
    const std::string& id() const;
    
    /** Returns 1, 0, or -1, if the operator is symmetric, nonsymmetric,
        or antisymmetric with respect to permutation of particles i and j */
    const int psymm(int i, int j) const;

  };
  
  //extern Operator TwoERep;
  //extern Operator TwoEDist;

  class IntegralType {

    // Operator
    Operator O_;
    
    // number of basis functions in bra for each particle >= 0
    vector<char> nbra_;
    // number of basis functions in ket for each particle >= 0
    vector<char> nket_;
    
    public:
    IntegralType(const Operator& O, const vector<char>& nbra, const vector<char>& nket);
    ~IntegralType();
    
  };


  class IntegralSetType {

    public:
    IntegralSetType();
    virtual ~IntegralSetType();
    
    /// Number of elements (>0 if static, 0 if determined at runtime)
    virtual int num_elems();

  };

  class SingleIntegral : public IntegralSetType {

    public:
    SingleIntegral();

  };

  class IntegralShell : public IntegralSetType {

    public:
    IntegralShell();

  };

  class IntegralShellVector : public IntegralSetType {

    public:
    IntegralShellVector();

  };

  class IntegralVectorShell : public IntegralSetType {

    public:
    IntegralVectorShell();

  };

  class RecurrenceRelation {

    public:
    RecurrenceRelation();
    ~RecurrenceRelation();
    
    virtual const std::string cpp_function_name() =0;
    virtual const std::string cpp_source_name() =0;
    virtual const std::string cpp_header_name() =0;
    virtual std::ostream& cpp_source(std::ostream&) =0;

  };

  class DataTypedRR : public RecurrenceRelation {

    IntegralSetType ST_;

    public:
    DataTypedRR(const IntegralSetType&);
    ~DataTypedRR();

  };


  template<typename Q> class QuantumSet {

    Q quanta_;

  public:
    QuantumSet(Q quanta);
    ~QuantumSet();

    /// Increment quantum number i
    void inc_quantum(unsigned int i) { quanta_.inc(i); };
    /// Decrement quantum number i
    void dec_quantum(unsigned int i) { quanta_.dec(i); };

  };

  template<typename T, unsigned int N> class QuantumNumbers {
    T qn_[N];

  public:
    QuantumNumbers(T qn[N]);
    ~QuantumNumbers();

    /// Increment quantum number i
    void inc(unsigned int i) { ++qn[i]; };
    /// Decrement quantum number i
    void dec(unsigned int i) {
      if (qn[i] == T(0))
        throw std::runtime_error("QuantumNumber::dec -- quantum number already zero");
      --qn[i];
    };
  };

  template<typename Q>
    QuantumSet<Q>::QuantumSet(Q quanta) :
    quanta_(quanta)
  {
  };

  template<typename Q>
    QuantumSet<Q>::~QuantumSet()
  {
  };

  template<typename T, unsigned int N>
    QuantumNumbers<T,N>::QuantumNumbers(T qn[N])
  {
      for(int i=0; i<N; i++) qn_[i] = qn[i];
  };

  template<typename T, unsigned int N>
    QuantumNumbers<T,N>::~QuantumNumbers()
  {
  };

  typedef QuantumNumbers<unsigned int, 1> ShellQuantumNumbers;
  typedef QuantumNumbers<unsigned int, 3> CartGaussQuantumNumbers;
  typedef QuantumSet< ShellQuantumNumbers > GaussShell;
  typedef QuantumSet< CartGaussQuantumNumbers > CartGauss;


  /// Set of basis functions
  class BFSet {

    public:
    BFSet() {};
    virtual ~BFSet() {};
    virtual unsigned int num_bf() const =0;

    /// Increments one of the quantum numbers
    virtual void inc() =0;
    /// Decrements one of the quantum numbers
    virtual void dec() =0;

  };

  /// Cartesian Gaussian Function
  class CGF : public BFSet {

    unsigned int qn_[3];

    public:
    CGF(unsigned int qn[3]);
    CGF(const CGF&);
    ~CGF();

    /// Returns the number of basis functions in the set (always 1)
    unsigned int num_bf() const { return 1; };

  };

  /// Cartesian Gaussian Shell
  class CGShell : public BFSet {

    unsigned int qn_[1];

    public:
    /// Default constructor creates an s-type shell
    CGShell();
    CGShell(unsigned int qn[1]);
    CGShell(const CGShell&);
    ~CGShell();
    CGShell& operator=(const CGShell&);

    /// Returns the number of basis functions in the set
    unsigned int num_bf() const { return (qn_[0]+1)*(qn_[0]+2)/2; };

    /// Implements purely virtual BFSet::dec
    void dec();

    /// Implements purely virtual BFSet::inc
    void inc();

    /// Print out the content
    void print(std::ostream& os = std::cout) const;

  };


  class DGVertex;
  /** Class DGArc describes arcs in a directed graph.
      Each arc connects vertex orig_ to vertex dest_. */
  class DGArc {

    const DGVertex* orig_;  // Where this Arc leavs
    const DGVertex* dest_;  // Where this Arc leads to

    public:
    DGArc(const DGVertex* orig_, const DGVertex* dest_);
    ~DGArc();
    
  };
  
  /** Class DGArcRel describes arcs in a directed graph which is
      represented by a relationship ArcRel. */
  // NOTE TO SELF (11/24/2004): need to implement checks on ArcRel
  // It obviously must implement some functions
  template <class ArcRel> class DGArcRel : public DGArc {

    ArcRel* rel_;     // Relationship described by the arc

    public:
    DGArcRel(const DGVertex* orig, const DGVertex* dest, const ArcRel* rel);
    ~DGArcRel();
    
  };

  template <class ArcRel>
    DGArcRel<ArcRel>::DGArcRel(const DGVertex* orig, const DGVertex* dest, const ArcRel* rel) :
    DGArc(orig,dest), rel_(rel)
    {
    };

  template <class ArcRel>
    DGArcRel<ArcRel>::~DGArcRel()
    {
    };

  /// This is a vertex of a Directed Graph (DG)
  class DGVertex {

    /// Arcs leaving this DGVertex
    vector<DGArc*> children_;
    /// We also need info about Arcs entering this DGVertex
    vector<DGArc*> parents_;

    public:
    DGVertex();
    DGVertex(const vector<DGArc*>& parents, const vector<DGArc*>& children);
    ~DGVertex();

  };

  template <unsigned int NP> class Oper {

    /// Described name
    const std::string descr_;
    /// short (<20 chars) ID label
    const std::string id_;

    public:
    Oper(const std::string& descr, const std::string& id);
    virtual ~Oper();

    /// Returns full description of the operator
    const std::string& descr() const;
    /// Returns short label for the operator
    const std::string& id() const;
    
    /** Returns 1, 0, or -1, if the operator is symmetric, nonsymmetric,
        or antisymmetric with respect to permutation of particles i and j */
    virtual const int psymm(int i, int j) const =0;

    static const unsigned int np = NP;

  };

  template <unsigned int NP>
    Oper<NP>::Oper(const std::string& descr, const std::string& id) :
    descr_(descr), id_(id)
    {
    }
  
  template <unsigned int NP>
    Oper<NP>::~Oper()
    {
    }
  
  template <unsigned int NP>
    const std::string&
    Oper<NP>::descr() const
    {
      return descr_;
    }
  
  template <unsigned int NP>
    const std::string&
    Oper<NP>::id() const
    {
      return id_;
    }

  class TwoERep : public Oper<2> {

    // symmetry W.R.T. permutation of each pair of particles
    // 1 -- symmetric, -1 -- antisymmetric, 0 -- nonsymmetric
    // stored as a lower triangle (diagonal not included)
    static const char psymm_[np*(np-1)/2];

    public:
    TwoERep();
    ~TwoERep();

    /** Returns 1, 0, or -1, if the operator is symmetric, nonsymmetric,
        or antisymmetric with respect to permutation of particles i and j */
    const int psymm(int i, int j) const;

  };

  /** This template is for an integral of 1 operator over (products of) one
      type of basis functions. Any Integral can be a DGVertex.
  */
  template <class Oper, class BFSet> class Integral : public DGVertex {

    static Oper O_;

    protected:
    vector<BFSet> bra_[Oper::np];
    vector<BFSet> ket_[Oper::np];
    
    public:
    Integral(const vector<BFSet> bra[Oper::np], const vector<BFSet> ket[Oper::np]);
    virtual ~Integral();

    /// Copy is permitted
    Integral& operator=(const Integral& source);

    /// Obtain BFsets
    const BFSet& bra(unsigned int particle, unsigned int i) const;
    const BFSet& ket(unsigned int particle, unsigned int i) const;

  };

  template <class Oper, class BFSet>
    Integral<Oper, BFSet>::Integral(const vector<BFSet> bra[Oper::np], const vector<BFSet> ket[Oper::np]) :
    DGVertex()
    {
      for(int p=0; p<Oper::np; p++) {
        bra_[p] = bra[p];
        ket_[p] = ket[p];
      }
    };

  template <class Oper, class BFSet>
    Integral<Oper, BFSet>::~Integral()
    {
    };

  template <class Oper, class BFSet>
    Integral<Oper, BFSet>&
    Integral<Oper, BFSet>::operator=(const Integral<Oper, BFSet>& source)
    {
      for(int p=0; p<Oper::np; p++) {
        bra_[p] = source.bra_[p];
        ket_[p] = source.ket_[p];
      }
    };


  template <class Oper, class BFSet>
    const BFSet&
    Integral<Oper, BFSet>::bra(unsigned int p, unsigned int i) const
    {
      return bra_[p][i];
    };

  template <class Oper, class BFSet>
    const BFSet&
    Integral<Oper, BFSet>::ket(unsigned int p, unsigned int i) const
    {
      return ket_[p][i];
    };


  /// Standard ERI shell quartet
  template <class BFSet> class TwoERep_2b2k : public Integral<TwoERep, BFSet> {

    unsigned int m_;  // auxiliary index

    public:
    TwoERep_2b2k(const vector<BFSet> bra[TwoERep::np], const vector<BFSet> ket[TwoERep::np], unsigned int m);

    unsigned int m() const { return m_; };

    void print(std::ostream& os = std::cout) const;
  };

  template <class BFSet>
    TwoERep_2b2k<BFSet>::TwoERep_2b2k(const vector<BFSet> bra[TwoERep::np], const vector<BFSet> ket[TwoERep::np], unsigned int m) :
    Integral<TwoERep, BFSet>(bra, ket), m_(m)
    {
      if (bra[0].size() != 1)
        throw std::runtime_error("TwoERep_2b2k<BFSet>::TwoERep_2b2k(bra[2],ket[2]) -- dimension of bra[0] must be 1");
      if (bra[1].size() != 1)
        throw std::runtime_error("TwoERep_2b2k<BFSet>::TwoERep_2b2k(bra[2],ket[2]) -- dimension of bra[1] must be 1");
      if (ket[0].size() != 1)
        throw std::runtime_error("TwoERep_2b2k<BFSet>::TwoERep_2b2k(bra[2],ket[2]) -- dimension of ket[0] must be 1");
      if (ket[1].size() != 1)
        throw std::runtime_error("TwoERep_2b2k<BFSet>::TwoERep_2b2k(bra[2],ket[2]) -- dimension of ket[1] must be 1");
    };

  template <class BFSet>
    void
    TwoERep_2b2k<BFSet>::print(std::ostream& os) const
    {
      os << "TwoERep_2b2k: m = " << m_ << endl;
      os << "shell bra1:" << endl;
      bra_[0][0].print(os);
      os << "shell bra2:" << endl;
      bra_[1][0].print(os);
      os << "shell ket1:" << endl;
      ket_[0][0].print(os);
      os << "shell ket2:" << endl;
      ket_[1][0].print(os);
    };

  /// VRR Recurrence Relation for ERI
  template <class BFSet> class VRR_ERI_2b2k : public RecurrenceRelation {

    const TwoERep_2b2k<BFSet>* target_;
    const TwoERep_2b2k<BFSet>* children_[5];

    public:
    VRR_ERI_2b2k(const TwoERep_2b2k<BFSet>*);
    ~VRR_ERI_2b2k();

    const std::string cpp_function_name() {};
    const std::string cpp_source_name() {};
    const std::string cpp_header_name() {};
    std::ostream& cpp_source(std::ostream&) {};

  };

  #include <vrr_eri_2b2k.h>

};

#endif

