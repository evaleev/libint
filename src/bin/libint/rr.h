

#ifndef _libint2_src_bin_libint_rr_h_
#define _libint2_src_bin_libint_rr_h_

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <exception.h>

using namespace std;


namespace libint2 {

  struct StaticDefinitions {
    static const unsigned int num_am_letters = 22;
    static const char am_letters[num_am_letters];
  };

  const int Libint2_DefaultVectorLength = 64;

  class RecurrenceRelation {

  public:
    RecurrenceRelation();
    virtual ~RecurrenceRelation();

    /** num_children() returns the actual number of children.
        For example, VRR for ERIs has 5 children on the right-hand side,
        however, for some ERI classes (ss|ps) the actual number may be
        smaller.
    */
    virtual const unsigned int num_children() const =0;
    
    virtual const std::string cpp_function_name() =0;
    virtual const std::string cpp_source_name() =0;
    virtual const std::string cpp_header_name() =0;
    virtual std::ostream& cpp_source(std::ostream&) =0;

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
    void inc(unsigned int i) { ++qn_[i]; };
    /// Decrement quantum number i
    void dec(unsigned int i) {
      if (qn_[i] == T(0))
        throw std::runtime_error("QuantumNumber::dec -- quantum number already zero");
      --qn_[i];
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
    virtual const std::string label() const =0;

    /// Increments one of the quantum numbers
    virtual void inc() =0;
    /// Decrements one of the quantum numbers
    virtual void dec() =0;

  };

  /// Cartesian Gaussian Function
  class CGF : public BFSet {

    unsigned int qn_[3];

  public:
    /// Default constructor makes an s-type Gaussian
    CGF();
    CGF(unsigned int qn[3]);
    CGF(const CGF&);
    CGF(const BFSet*);
    ~CGF();

    /// As far as SetIterator is concerned, CGF is a set of one CGF
    typedef CGF iter_type;

    /// Return a compact label
    const std::string label() const;
    /// Returns the number of basis functions in the set (always 1)
    unsigned int num_bf() const { return 1; };

    /// Returns the angular momentum
    unsigned int qn(unsigned int xyz) const;

    /// Comparison operator
    bool operator==(const CGF&) const;
    
    /// Decrement one of quantum numbers
    void dec();
    /// Increment one of quantum numbers
    void inc();

    /// Print out the content
    void print(std::ostream& os = std::cout) const;
    
  };

  /// Cartesian Gaussian Shell
  class CGShell : public BFSet {

    unsigned int qn_[1];

  public:
    /// Default constructor creates an s-type shell
    CGShell();
    CGShell(unsigned int qn[1]);
    CGShell(const CGShell&);
    CGShell(const BFSet*);
    ~CGShell();
    CGShell& operator=(const CGShell&);

    /// As far as SetIterator is concerned, CGShell is a set of one CGF
    typedef CGF iter_type;

    /// Return a compact label
    const std::string label() const;
    /// Returns the number of basis functions in the set
    unsigned int num_bf() const { return (qn_[0]+1)*(qn_[0]+2)/2; };

    /// Returns the angular momentum
    unsigned int qn() const { return qn_[0]; }

    /// Comparison operator
    bool operator==(const CGShell&) const;

    /// Implements purely virtual BFSet::dec, may throw InvalidDecrement
    void dec();
    /// Implements purely virtual BFSet::inc
    void inc() throw();

    /// Print out the content
    void print(std::ostream& os = std::cout) const;

  };


  class DGVertex;
  /** Class DGArc describes arcs in a directed graph.
      Each arc connects vertex orig_ to vertex dest_. */
  class DGArc {

    DGVertex* orig_;  // Where this Arc leavs
    DGVertex* dest_;  // Where this Arc leads to

  public:
    DGArc(DGVertex* orig_, DGVertex* dest_);
    ~DGArc();

    DGVertex* orig() const { return orig_; };
    DGVertex* dest() const { return dest_; };
    
  };
  
  /** Class DGArcRel describes arcs in a directed graph which is
      represented by a relationship ArcRel. */
  // NOTE TO SELF (11/24/2004): need to implement checks on ArcRel
  // It obviously must implement some functions
  template <class ArcRel> class DGArcRel : public DGArc {

    ArcRel* rel_;     // Relationship described by the arc

  public:
    DGArcRel(DGVertex* orig, DGVertex* dest, ArcRel* rel);
    ~DGArcRel();
    
  };

  template <class ArcRel>
    DGArcRel<ArcRel>::DGArcRel(DGVertex* orig, DGVertex* dest, ArcRel* rel) :
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

    // Whether this is a "target" vertex, i.e. the target of a calculation
    bool target_;
    // If set to true -- traversal has started and add_entry... cannot be called
    bool can_add_arcs_;
    /// add_entry_arc(arc) adds arc as an arc connecting parents to this vertex
    void add_entry_arc(DGArc*);
    /// del_entry_arc(arc) removes arc as an arc connecting parents to this vertex
    void del_entry_arc(DGArc*);

    ////////
    // These members used in traversal algorithms
    ////////

    // num_tagged_arcs keeps track of how many entry arcs have been tagged during traversal
    unsigned int num_tagged_arcs_;
    /// Which DGVertex to be computed before this vertex (0, if this is the first vertex)
    DGVertex* precalc_;
    /// Which DGVertex to be computed after this vertex (0, if this is the last vertex)
    DGVertex* postcalc_;

  public:
    DGVertex();
    DGVertex(const vector<DGArc*>& parents, const vector<DGArc*>& children);
    virtual ~DGVertex();

    /// make_a_target() marks this vertex as a target
    void make_a_target();
    /// is_a_target() returns true if this vertex is a target
    const bool is_a_target() const { return target_;};
    /** add_exit_arc(arc) adds arc as an arc connecting to children of this vertex.
        Thus, arcs are owned by their PARENTS.
      */
    void add_exit_arc(DGArc*);
    /// returns the number of parents
    const unsigned int num_entry_arcs() const;
    /// returns ptr to i-th parent
    DGArc* entry_arc(unsigned int) const;
    /// returns the number of children
    const unsigned int num_exit_arcs() const;
    /// returns ptr to i-th child
    DGArc* exit_arc(unsigned int) const;

    /** apply_rr() applies the optimal recurrence relation to this particular DGVertex.
        The concrete class must implement this.
    */
    //virtual RecurrenceRelation* apply_rr() =0;

    /** equiv(const DGVertex* aVertex) returns true if this vertex is
        equivalent to *aVertex.
    */
    virtual bool equiv(const DGVertex*) const =0;

    /** print(std::ostream&) prints out comment-style info vertex
    */
    virtual void print(std::ostream& os = std::cout) const =0;

    /// prepare_to_traverse() must be called before traversal of the graph starts
    void prepare_to_traverse();
    /// tag() tags the vertex and returns the total number of tags this vertex has received
    const unsigned int tag();
    /// Returns pointer to vertex to be computed before this vertex, 0 if this is the first vertex
    DGVertex* precalc() const { return precalc_; };
    /// Returns pointer to vertex to be computed after this vertex, 0 if this is the last vertex
    DGVertex* postcalc() const { return postcalc_; };
    /// Sets precalc
    void set_precalc(DGVertex* precalc) { precalc_ = precalc; };
    /// Sets postcalc
    void set_postcalc(DGVertex* postcalc) { postcalc_ = postcalc; };

    /// Resets the vertex, releasing all arcs
    void reset();
  };

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
    */
  class OperSet {
    public:
      virtual ~OperSet() {};

      /// Returns full description of the operator
      virtual const std::string& descr() const =0;
      /// Returns short label for the operator
      virtual const std::string& id() const =0;

      /** Returns 1, 0, or -1, if each operator in the set is symmetric, nonsymmetric,
        or antisymmetric with respect to permutation of particles i and j */
      virtual const int psymm(int i, int j) const =0;

      /// Number of operators in the set
      virtual const unsigned int num_oper() const =0;
    };
  

  template <class Props>
    class Oper : public OperSet {
    public:
      typedef Props Properties;      
      virtual ~Oper();

      /// Returns full description of the operator
      const std::string& descr() const;
      /// Returns short label for the operator
      const std::string& id() const;

    protected:
      /// The only declared constructor is only useable by derived classes
      Oper(const std::string& descr, const std::string& id);

    private:
      /// Described name
      const std::string descr_;
      /// short (<20 chars) ID label
      const std::string id_;
  };

  template <class Props>
    Oper<Props>::Oper(const std::string& descr, const std::string& id) :
    OperSet(), descr_(descr), id_(id)
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
    Oper<Props>::id() const
    {
      return id_;
    }

  typedef OperatorProperties<2,true> TwoPRep_Props;
  class TwoERep : public Oper<TwoPRep_Props> {
  public:
    /// TwoERep is not a set
    typedef TwoERep iter_type;
    const unsigned int num_oper() const { return 1; };
  
    TwoERep();
    TwoERep(const OperSet*);
    ~TwoERep();

    /** Returns 1, 0, or -1, if the operator is symmetric, nonsymmetric,
        or antisymmetric with respect to permutation of particles i and j */
    const int psymm(int i, int j) const;

  private:
    // symmetry W.R.T. permutation of each pair of particles
    // 1 -- symmetric, -1 -- antisymmetric, 0 -- nonsymmetric
    // stored as a lower triangle (diagonal not included)
    static const char psymm_[Properties::np*(Properties::np-1)/2];

  };

  /** This template is for an integral of 1 operator over (products of) one
      type of basis functions. No instances of Integral can be created --
      this class is intended as a base class.
  */
  template <class Oper, class BFSet> class Integral {

    static Oper O_;

    //
    // All integrals are Singletons by nature, therefore they must be treated as such
    // 1) No public constructors are provided
    // 2) protected members are provided to implement Singleton-type functionality
    //
    Integral();
    Integral(const Integral&);
    // Copy is not permitted
    Integral& operator=(const Integral& source);

  protected:
    // Basic Integral constructor
    Integral(const vector<BFSet> bra[Oper::Properties::np], const vector<BFSet> ket[Oper::Properties::np]);

    vector<BFSet> bra_[Oper::Properties::np];
    vector<BFSet> ket_[Oper::Properties::np];
    
  public:
    /** No constructors are public since this is a singleton-like quantity.
        Instead, access is provided through derived class's Instance().
     */
    virtual ~Integral();

    /// Equivalence operator
    virtual bool equiv(const Integral*) const;

    /// Obtain BFsets
    const BFSet& bra(unsigned int particle, unsigned int i) const;
    const BFSet& ket(unsigned int particle, unsigned int i) const;

  };

  template <class Oper, class BFSet>
    Integral<Oper, BFSet>::Integral(const vector<BFSet> bra[Oper::Properties::np], const vector<BFSet> ket[Oper::Properties::np])
    {
      for(int p=0; p<Oper::Properties::np; p++) {
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
      for(int p=0; p<Oper::Properties::np; p++) {
        bra_[p] = source.bra_[p];
        ket_[p] = source.ket_[p];
      }
    };

  template <class Oper, class BFSet>
    bool
    Integral<Oper, BFSet>::equiv(const Integral<Oper,BFSet>* a) const
    {
      bool equiv = true;
      for(int p=0; p<Oper::Properties::np; p++) {
        equiv = (equiv && bra_[p] == a->bra_[p]);
        equiv = (equiv && ket_[p] == a->ket_[p]);
      }
      //if (equiv)
      //  cout << "Integral<Oper, BFSet>::equiv() returned true" << endl;
      return equiv;
    }

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
  template <class BFSet> class TwoERep_2b2k : public Integral<TwoERep, BFSet>, public DGVertex {

    unsigned int m_;  // auxiliary index

    // Default and copy constructors are not allowed
    TwoERep_2b2k();
    TwoERep_2b2k(const TwoERep_2b2k&);

    // This constructor is private since all Integral's are Singletons. Use Instance instead.
    TwoERep_2b2k(const vector<BFSet> bra[TwoERep::Properties::np], const vector<BFSet> ket[TwoERep::Properties::np], unsigned int m);
    // stack_ of pointers to objects used to check whether an object already exists
    static vector< TwoERep_2b2k* > stack_;

  public:
    /// Returns a pointer to a unique instance, a la Singleton
    static TwoERep_2b2k* Instance(const vector<BFSet> bra[TwoERep::Properties::np], const vector<BFSet> ket[TwoERep::Properties::np], unsigned int m);
    
    unsigned int m() const { return m_; };

    /// Overload of DGVertex's equiv
    bool equiv(const DGVertex*) const;

    void print(std::ostream& os = std::cout) const;
  };

  template <class BFSet>
    vector< TwoERep_2b2k<BFSet>* > TwoERep_2b2k<BFSet>::stack_(0);

  template <class BFSet>
    TwoERep_2b2k<BFSet>::TwoERep_2b2k(const vector<BFSet> bra[TwoERep::Properties::np], const vector<BFSet> ket[TwoERep::Properties::np], unsigned int m) :
    Integral<TwoERep, BFSet>(bra, ket), DGVertex(), m_(m)
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
    TwoERep_2b2k<BFSet>*
    TwoERep_2b2k<BFSet>::Instance(const vector<BFSet> bra[TwoERep::Properties::np], const vector<BFSet> ket[TwoERep::Properties::np], unsigned int m)
    {
      TwoERep_2b2k* const this_int = new TwoERep_2b2k<BFSet>(bra,ket,m);
      int stack_size = stack_.size();
      for(int i=0; i<stack_size; i++) {
        if (this_int->equiv(stack_[i])) {
          delete this_int;
          return stack_[i];
        }
      }
      stack_.push_back(this_int);
      return this_int;
    };

  template <class BFSet>
    bool
    TwoERep_2b2k<BFSet>::equiv(const DGVertex* a) const
    {
      // check the type first
      const TwoERep_2b2k<BFSet>* a_cast = dynamic_cast< const TwoERep_2b2k<BFSet>* >(a);
      if (!a_cast)
        return false;

      bool result = Integral<TwoERep, BFSet>::equiv(a_cast) && (m_ == a_cast->m_);
      return result;
    }

  template <class BFSet>
    void
    TwoERep_2b2k<BFSet>::print(std::ostream& os) const
    {
      os << "TwoERep_2b2k: (" << Integral<TwoERep,BFSet>::bra(0,0).label() << " " << Integral<TwoERep,BFSet>::ket(0,0).label()
         << " | " << Integral<TwoERep,BFSet>::bra(1,0).label() << " " << Integral<TwoERep,BFSet>::ket(1,0).label() << ")^{" << m_ <<"}" << endl;
    };

  /** VRR Recurrence Relation for 2-e ERI. part specifies for which particle
      the angular momentum is raised. bool bra specifies whether the angular momentum
      is raised in bra (true) or ket (false). Class ERI specifies which particular implementation
      of ERI to use.
   */
  template <template <class> class ERI, class BFSet, int part, bool bra> class VRR_ERI_2b2k : public RecurrenceRelation {

    static const unsigned int nchild_ = 5;

    ERI<BFSet>* target_;
    ERI<BFSet>* children_[nchild_];

    unsigned int num_actual_children_;

  public:
    VRR_ERI_2b2k(ERI<BFSet>*);
    ~VRR_ERI_2b2k();

    typedef ERI<BFSet> TargetType;

    const unsigned int num_children() const { return num_actual_children_; };
    /// target() returns points to the i-th child
    ERI<BFSet>* target() { return target_; };
    /// child(i) returns points i-th child
    ERI<BFSet>* child(unsigned int i);

    const std::string cpp_function_name() {};
    const std::string cpp_source_name() {};
    const std::string cpp_header_name() {};
    std::ostream& cpp_source(std::ostream&) {};

  };

  /** A generic HRR Recurrence Relation. Int is the integral class. part specifies for which particle
      the angular momentum is shifted.
   */
  template <template <class> class I, class BFSet, int part> class HRR_ERI_2b2k : public RecurrenceRelation {

    static const unsigned int nchild_ = 2;

    I<BFSet>* target_;
    I<BFSet>* children_[nchild_];

    unsigned int num_actual_children_;

  public:
    HRR_ERI_2b2k(I<BFSet>*);
    ~HRR_ERI_2b2k();

    typedef I<BFSet> TargetType;

    const unsigned int num_children() const { return num_actual_children_; };
    /// target() returns points to the i-th child
    I<BFSet>* target() { return target_; };
    /// child(i) returns points i-th child
    I<BFSet>* child(unsigned int i);

    const std::string cpp_function_name() {};
    const std::string cpp_source_name() {};
    const std::string cpp_header_name() {};
    std::ostream& cpp_source(std::ostream&) {};

  };

  typedef enum {
    InBra=0, InKet=1
  } FunctionPosition;
  typedef enum {
    BraToKet=0, KetToBra=1
  } FunctionMovement;
  
  /** VRR Recurrence Relation for 2-e ERI. part specifies for which particle
      the angular momentum is raised. bool bra specifies whether the angular momentum
      is raised in bra (true) or ket (false). Class ERI specifies which particular implementation
      of ERI to use.
   */
  template <template <class> class ERI, class BFSet, int part, FunctionPosition where>
    class VRR_11_TwoPRep_11 : public RecurrenceRelation {

    static const unsigned int nchild_ = 5;

    ERI<BFSet>* target_;
    ERI<BFSet>* children_[nchild_];

    unsigned int num_actual_children_;

  public:
    VRR_11_TwoPRep_11(ERI<BFSet>*);
    ~VRR_11_TwoPRep_11();

    typedef ERI<BFSet> TargetType;

    const unsigned int num_children() const { return num_actual_children_; };
    /// target() returns points to the i-th child
    ERI<BFSet>* target() { return target_; };
    /// child(i) returns points i-th child
    ERI<BFSet>* child(unsigned int i);

    const std::string cpp_function_name() {};
    const std::string cpp_source_name() {};
    const std::string cpp_header_name() {};
    std::ostream& cpp_source(std::ostream&) {};

  };

  /** A generic Horizontal Recurrence Relation:

      |a b) = |a+1 b-1) + AB |a b-1)

      Int is the integral class. part specifies for which particle
      the angular momentum is shifted. Function a is assumed to gain quanta,
      function a gains quanta. loc_a and loc_b specify where
      functions a and b are located (bra or ket). pos_a and pos_b
      specify which function to be used (usually pos_a and pos_b are set
      to 0 to refer to the first function for this particle in this location).
      
   */
  template <template <class> class I, class BFSet, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    class HRR : public RecurrenceRelation {

    static const unsigned int nchild_ = 2;

    I<BFSet>* target_;
    I<BFSet>* children_[nchild_];

    unsigned int num_actual_children_;

    void oper_checks() const;

  public:
    HRR(I<BFSet>*);
    ~HRR();

    typedef I<BFSet> TargetType;

    const unsigned int num_children() const { return num_actual_children_; };
    /// target() returns points to the i-th child
    I<BFSet>* target() { return target_; };
    /// child(i) returns points i-th child
    I<BFSet>* child(unsigned int i);

    const std::string cpp_function_name() {};
    const std::string cpp_source_name() {};
    const std::string cpp_header_name() {};
    std::ostream& cpp_source(std::ostream&) {};

  };


};

#include <vrr_eri_2b2k.h>
#include <hrr_eri_2b2k.h>
#include <vrr_11_twoprep_11.h>
#include <hrr.h>
//#include <shell_to_ints.h>
//#include <iter.h>

#endif

