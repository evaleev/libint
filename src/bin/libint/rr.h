
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <exception.h>
#include <bfset.h>
#include <smart_ptr.h>
#include <polyconstr.h>
#include <singl_stack.h>

#ifndef _libint2_src_bin_libint_rr_h_
#define _libint2_src_bin_libint_rr_h_

using namespace std;

namespace libint2 {

  class DGVertex;
  class CodeContext;

  /**
     RecurrenceRelation describes all recurrence relations
  */
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
    /**
       Returns true is the type of target and all children are exactly the same
    */
    virtual bool invariant_type() const =0;

    /**
      label() returns a unique, short, descriptive label of this RR
      (e.g. "VRR A (p s | 1/r_{12} | d s )" for Obara-Saika recurrence relation
      applied to center A to compute (ps|ds) ERI)
    */
    virtual std::string label() const =0;
    
    /// Produce an expression which implements this RR by calling the specialized function
    virtual void spfunction_call(const SafePtr<CodeContext>& context, std::ostream& os) const =0;
    
    /// Generate code for the recurrence relation
    void generate_code(const SafePtr<CodeContext>& context,
                       std::ostream& decl, std::ostream& def);

    /// Return the number of FLOPs per this recurrence relation
    virtual unsigned int nflops() const { return 0; }
    
    virtual const std::string cpp_function_name() =0;
    virtual const std::string cpp_source_name() =0;
    virtual const std::string cpp_header_name() =0;
    virtual std::ostream& cpp_source(std::ostream&) =0;

  };


  /** RRTactic describes an object that specifies a tactic of how to apply
    recurrence relation
    */

  typedef enum {
    InBra=0, InKet=1
  } FunctionPosition;
  typedef enum {
    BraToKet=0, KetToBra=1
  } FunctionMovement;
  
  /** RRStack implements a stack of RecurrenceRelation's which can only hold
      one instance of a given RR. RecurrenceRelation::label() is used for hashing
    */
  class RRStack : public SingletonStack<RecurrenceRelation,std::string>
  {
    public:
    typedef SingletonStack<RecurrenceRelation,std::string> parent_type;
    typedef parent_type::data_type data_type;
    typedef parent_type::iter_type iter_type;
    typedef parent_type::citer_type citer_type;
    
    RRStack() : parent_type(&RecurrenceRelation::label) {}
    ~RRStack() {}
    
    /// adds content of rrs to this stack
    void add(const SafePtr<RRStack>& rrs);
  };
  
};

//#include <vrr_11_twoprep_11.h>
//#include <hrr.h>
//#include <shell_to_ints.h>
//#include <iter.h>

#endif

