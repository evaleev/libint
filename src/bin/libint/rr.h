
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <exception.h>
#include <smart_ptr.h>
#include <polyconstr.h>

#ifndef _libint2_src_bin_libint_rr_h_
#define _libint2_src_bin_libint_rr_h_

using namespace std;


namespace libint2 {

  struct StaticDefinitions {
    static const unsigned int num_am_letters = 22;
    static const char am_letters[num_am_letters];
  };

  const int Libint2_DefaultVectorLength = 64;

  class DGVertex;
  
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
    /// Returns i-th child
    virtual SafePtr<DGVertex> rr_child(unsigned int i) const =0;
    /// Returns the target
    virtual SafePtr<DGVertex> rr_target() const =0;
    /** num_expr() returns the actual number of expressions.
        For example, VRR for ERIs has up to 6 expressions on the right-hand side,
    */
    virtual const unsigned int num_expr() const =0;
    /// Returns i-th expression
    virtual SafePtr<DGVertex> rr_expr(unsigned int i) const =0;
    /**
       Returns true is this recurrence relation is simple enough to optimize away.
       As a result of such optimization, standalone function will NOT be 
       generated for this recurrence relation. Instead, it's source will be
       inlined and optimized.
    */
    virtual bool is_simple() const =0;

    /// Return the number of FLOPs per this recurrence relation
    virtual unsigned int nflops() const { return 0; }

    virtual const std::string cpp_function_name() =0;
    virtual const std::string cpp_source_name() =0;
    virtual const std::string cpp_header_name() =0;
    virtual std::ostream& cpp_source(std::ostream&) =0;

  };

  /** Set of basis functions. Sets must be constructable using
      SafePtr<BFSet> or SafePtr<ConstructablePolymorphically>.
  */
  class BFSet : public ConstructablePolymorphically {

  public:
    virtual ~BFSet() {}
    virtual unsigned int num_bf() const =0;
    virtual const std::string label() const =0;

  protected:
    BFSet() {}

  };

  /** Set of basis functions with incrementable/decrementable quantum numbers.
      Sets must be constructable using SafePtr<BFSet> or SafePtr<ConstructablePolymorphically>.
  */
  class IncableBFSet : public BFSet {

  public:
    virtual ~IncableBFSet() {}

    /// Increment i-th quantum number. Do nothing if i is outside the allowed range
    virtual void inc(unsigned int i) throw() =0;
    /// Decrements i-th quantum number. Do nothing is i is outside the allowed range
    virtual void dec(unsigned int i) =0;
    /// Returns true if all quanta are 0
    virtual bool zero() const =0;
    
  protected:
    IncableBFSet() {}

  };

  /// Cartesian Gaussian Function
  class CGF : public IncableBFSet {

    unsigned int qn_[3];

  public:
    /// As far as SetIterator is concerned, CGF is a set of one CGF
    typedef CGF iter_type;
    typedef IncableBFSet parent_type;

    /// Default constructor makes an s-type Gaussian
    CGF();
    CGF(unsigned int qn[3]);
    CGF(const CGF&);
    CGF(const SafePtr<CGF>&);
    CGF(const SafePtr<parent_type>&);
    CGF(const SafePtr<ConstructablePolymorphically>&);
    ~CGF();

    /// Return a compact label
    const std::string label() const;
    /// Returns the number of basis functions in the set (always 1)
    unsigned int num_bf() const { return 1; };

    /// Returns the angular momentum
    unsigned int qn(unsigned int xyz) const;

    /// Comparison operator
    bool operator==(const CGF&) const;
    
    /// Implements purely virtual IncableBFSet::dec, may throw InvalidDecrement
    void dec(unsigned int i);
    /// Implements purely virtual IncableBFSet::inc
    void inc(unsigned int i) throw();
    /// Implements IncableBFSet::zero()
    bool zero() const;

    /// Print out the content
    void print(std::ostream& os = std::cout) const;
    
  };

  /// Cartesian Gaussian Shell
  class CGShell : public IncableBFSet {

    unsigned int qn_[1];

  public:
    /// As far as SetIterator is concerned, CGShell is a set of one CGF
    typedef CGF iter_type;
    typedef IncableBFSet parent_type;

    /// Default constructor creates an s-type shell
    CGShell();
    CGShell(unsigned int qn[1]);
    CGShell(const CGShell&);
    CGShell(const SafePtr<CGShell>&);
    CGShell(const SafePtr<parent_type>&);
    CGShell(const SafePtr<ConstructablePolymorphically>&);
    ~CGShell();
    CGShell& operator=(const CGShell&);

    /// Return a compact label
    const std::string label() const;
    /// Returns the number of basis functions in the set
    unsigned int num_bf() const { return (qn_[0]+1)*(qn_[0]+2)/2; };

    /// Returns the angular momentum
    unsigned int qn(unsigned int xyz=0) const { return qn_[0]; }

    /// Comparison operator
    bool operator==(const CGShell&) const;

    /// Implements purely virtual IncableBFSet::dec, may throw InvalidDecrement
    void dec(unsigned int i);
    /// Implements purely virtual IncableBFSet::inc
    void inc(unsigned int i) throw();
    /// Implements IncableBFSet::zero()
    bool zero() const;

    /// Print out the content
    void print(std::ostream& os = std::cout) const;

  };

  /**
     TrivialBFSet<T> defines static member result, which is true if T
     is a basis function set consisting of 1 function
  */
  template <class T>
    struct TrivialBFSet;
  template <>
    struct TrivialBFSet<CGShell> {
      static const bool result = false;
    };
  template <>
    struct TrivialBFSet<CGF> {
      static const bool result = true;
    };


  /** RRTactic describes an object that specifies a tactic of how to apply
    recurrence relation
    */

  class DGVertex;
  /** Class DGArc describes arcs in a directed graph.
      Each arc connects vertex orig_ to vertex dest_. */
  class DGArc {

    SafePtr<DGVertex> orig_;  // Where this Arc leavs
    SafePtr<DGVertex> dest_;  // Where this Arc leads to

  public:
    DGArc(const SafePtr<DGVertex>& orig, const SafePtr<DGVertex>& dest);
    ~DGArc() {}

    SafePtr<DGVertex> orig() const { return orig_; }
    SafePtr<DGVertex> dest() const { return dest_; }

    /// Print out the arc
    virtual void print(std::ostream&) const =0;

  };
  
  /** Class DGArcDirect describes arcs that does not correspond to any relationship.
      Each arc connects vertex orig_ to vertex dest_. */
  class DGArcDirect : public DGArc {

  public:
    DGArcDirect(const SafePtr<DGVertex>& orig, const SafePtr<DGVertex>& dest) : DGArc(orig,dest) {}
    ~DGArcDirect() {}

    /// Overload of DGArc::print()
    void print(std::ostream& os) const
      {
        os << "DGArcDirect: connects " << orig().get() << " to " << dest().get();
      }
  };
  
  /** Class DGArcRR describes arcs correspond to recurrence relations.
      Each arc connects vertex orig_ to vertex dest_. */
  class DGArcRR : public DGArc {

  public:
    ~DGArcRR() {}

    /// rr() returns pointer to the RecurrenceRelation describing the arc
    virtual SafePtr<RecurrenceRelation> rr() const =0;

  protected:
    DGArcRR(const SafePtr<DGVertex>& orig, const SafePtr<DGVertex>& dest);

  };
  
  /** Class DGArcRel describes arcs in a directed graph which is
      represented by a relationship ArcRel. */
  // NOTE TO SELF (11/24/2004): need to implement checks on ArcRel
  // It obviously must implement some functions
  template <class ArcRel> class DGArcRel : public DGArcRR {

    SafePtr<ArcRel> rel_;     // Relationship described by the arc

  public:
    DGArcRel(const SafePtr<DGVertex>& orig, const SafePtr<DGVertex>& dest,
	     const SafePtr<ArcRel>& rel);
    ~DGArcRel();

    /// Implementation of DGArcRR::rr()
    SafePtr<RecurrenceRelation> rr() const { return dynamic_pointer_cast<RecurrenceRelation,ArcRel>(rel_); }
    /// Overload of DGArc::print()
    void print(std::ostream& os) const
      {
        os << "DGArcRel<T>: connects " << orig().get() << " to " << dest().get();
      }
    
  };

  template <class ArcRel>
    DGArcRel<ArcRel>::DGArcRel(const SafePtr<DGVertex>& orig, const SafePtr<DGVertex>& dest,
			       const SafePtr<ArcRel>& rel) :
    DGArcRR(orig,dest), rel_(rel)
    {
    };

  template <class ArcRel>
    DGArcRel<ArcRel>::~DGArcRel()
    {
    };

  /// This is a vertex of a Directed Graph (DG)
  class DGVertex {

    /// vertex label
    std::string label_;
    
    /// Arcs leaving this DGVertex
    vector< SafePtr<DGArc> > children_;
    /// We also need info about Arcs entering this DGVertex
    vector< SafePtr<DGArc> > parents_;

    // Whether this is a "target" vertex, i.e. the target of a calculation
    bool target_;
    // If set to true -- traversal has started and add_entry... cannot be called
    bool can_add_arcs_;
    /// add_entry_arc(arc) adds arc as an arc connecting parents to this vertex
    void add_entry_arc(const SafePtr<DGArc>&);
    /// del_entry_arc(arc) removes arc as an arc connecting parents to this vertex
    void del_entry_arc(const SafePtr<DGArc>&);

    ////////
    // These members used in traversal algorithms
    ////////

    // num_tagged_arcs keeps track of how many entry arcs have been tagged during traversal
    unsigned int num_tagged_arcs_;
    /// Which DGVertex to be computed before this vertex (0, if this is the first vertex)
    SafePtr<DGVertex> precalc_;
    /// Which DGVertex to be computed after this vertex (0, if this is the last vertex)
    SafePtr<DGVertex> postcalc_;

  public:
    DGVertex();
    DGVertex(const vector<SafePtr<DGArc> >& parents, const vector<SafePtr<DGArc> >& children);
    virtual ~DGVertex();

    /// make_a_target() marks this vertex as a target
    void make_a_target();
    /// is_a_target() returns true if this vertex is a target
    const bool is_a_target() const { return target_;};
    /** add_exit_arc(arc) adds arc as an arc connecting to children of this vertex.
        Thus, arcs are owned by their PARENTS.
      */
    void add_exit_arc(const SafePtr<DGArc>&);
    /** del_exit_arc(arc) removes arc c (from this and corresponding child)
      */
    void del_exit_arc(const SafePtr<DGArc>&);
    /// returns the number of parents
    const unsigned int num_entry_arcs() const;
    /// returns ptr to i-th parent
    SafePtr<DGArc> entry_arc(unsigned int) const;
    /// returns the number of children
    const unsigned int num_exit_arcs() const;
    /// returns ptr to i-th child
    SafePtr<DGArc> exit_arc(unsigned int) const;

    /** apply_rr() applies the optimal recurrence relation to this particular DGVertex.
        The concrete class must implement this.
    */
    //virtual RecurrenceRelation* apply_rr() =0;

    /** equiv(const DGVertex* aVertex) returns true if this vertex is
        equivalent to *aVertex.
    */
    virtual bool equiv(const SafePtr<DGVertex>&) const =0;
    
    /** precomputed() returns whether this DGVertex is precomputed
    */
    virtual bool precomputed() const =0;

    /** Returns the amount of memory (in floating-point words) to be allocated for the vertex.
      */
    virtual const unsigned int size() const =0;
    
    /** print(std::ostream&) prints out comment-style info vertex
    */
    virtual void print(std::ostream& os = std::cout) const =0;
    /// returns the label
    const std::string& label() const { return label_;}
    /// sets the label
    void set_label(const std::string& label);

    /// prepare_to_traverse() must be called before traversal of the graph starts
    void prepare_to_traverse();
    /// tag() tags the vertex and returns the total number of tags this vertex has received
    const unsigned int tag();
    /// Returns pointer to vertex to be computed before this vertex, 0 if this is the first vertex
    SafePtr<DGVertex> precalc() const { return precalc_; };
    /// Returns pointer to vertex to be computed after this vertex, 0 if this is the last vertex
    SafePtr<DGVertex> postcalc() const { return postcalc_; };
    /// Sets precalc
    void set_precalc(const SafePtr<DGVertex>& precalc) { precalc_ = precalc; };
    /// Sets postcalc
    void set_postcalc(const SafePtr<DGVertex>& postcalc) { postcalc_ = postcalc; };

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
      virtual const int psymm(int i, int j) const =0;

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
    
  typedef OperatorProperties<2,true> TwoPRep_Props;
  class TwoERep : public Oper<TwoPRep_Props> {
  public:
    typedef Oper<TwoPRep_Props> parent_type;
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
    const int psymm(int i, int j) const;

  private:
    // symmetry W.R.T. permutation of each pair of particles
    // 1 -- symmetric, -1 -- antisymmetric, 0 -- nonsymmetric
    // stored as a lower triangle (diagonal not included)
    static const char psymm_[Properties::np*(Properties::np-1)/2];

  };


  typedef enum {
    InBra=0, InKet=1
  } FunctionPosition;
  typedef enum {
    BraToKet=0, KetToBra=1
  } FunctionMovement;
  
};

#include <vrr_11_twoprep_11.h>
#include <hrr.h>
//#include <shell_to_ints.h>
//#include <iter.h>

#endif

