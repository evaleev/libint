
#ifndef _libint2_src_bin_libint_rr_h_
#define _libint2_src_bin_libint_rr_h_

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
#include <code.h>
#include <default_params.h>
#include <util_types.h>

using namespace std;

namespace libint2 {

  class DGVertex;
  class CodeContext;
  class ImplicitDimensions;
  class DirectedGraph;
  template <typename V> class AlgebraicOperator;

  /**
     RecurrenceRelation describes all recurrence relations
  */
  class RecurrenceRelation : public EnableSafePtrFromThis<RecurrenceRelation> {
  public:
    typedef RecurrenceRelation this_type;

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

    /** Numerical expression of a recurrence relation is always expressed as
        an AlgebraicOperator<DGVertex> */
    typedef AlgebraicOperator<DGVertex> ExprType;
    /// Returns the expression
    const SafePtr<ExprType>& rr_expr() const { return expr_; }

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
    virtual bool invariant_type() const;

    /**
      label() returns a unique, short, descriptive label of this RR
      (e.g. "VRR A (p s | 1/r_{12} | d s )" for Obara-Saika recurrence relation
      applied to center A to compute (ps|ds) ERI)
    */
    const std::string& label() const { if (label_.empty()) label_ = generate_label(); return label_; }
    
    /**
      description() returns a verbose description of this RR
    */
    virtual const std::string& description() const;
    
    /// Generate declaration and definition for the recurrence relation
    virtual void generate_code(const SafePtr<CodeContext>& context,
                               const SafePtr<ImplicitDimensions>& dims,
                               const std::string& funcname,
                               std::ostream& decl, std::ostream& def);

    /// Generate declaration and definition for the recurrence relation
    /// using generic code (typically, a manually written code)
    virtual void generate_generic_code(const SafePtr<CodeContext>& context,
                                       const SafePtr<ImplicitDimensions>& dims,
                                       const std::string& funcname,
                                       std::ostream& decl, std::ostream& def);

    /// Generate a callback for this recurrence relation
    virtual std::string spfunction_call(const SafePtr<CodeContext>& context,
                                        const SafePtr<ImplicitDimensions>& dims) const;

    /// Return the number of FLOPs per this recurrence relation
    unsigned int nflops() const { return nflops_; }
    
    virtual const std::string cpp_function_name() =0;
    virtual const std::string cpp_source_name() =0;
    virtual const std::string cpp_header_name() =0;
    virtual std::ostream& cpp_source(std::ostream&) =0;

    /// RecurrenceRelation is managed by SingletonStack but doesn't need to keep track of instance ID
    void inst_id(const SingletonStack<RecurrenceRelation,string>::InstanceID& i) {}

    protected:
    unsigned int nflops_;
    mutable std::string label_;
    SafePtr<ExprType> expr_;
    /// Adds a (or -a, if minus = -1) to expr_.
    void add_expr(const SafePtr<ExprType>& a, int minus=1);
    /// Generates the label
    virtual std::string generate_label() const =0;
    /// Registers with the stack
    template <class RR> bool register_with_rrstack() const;

    private:

    /** used by generate_code to create a (new) computation graph that computes sets of integrals using the RR
    */
    SafePtr<DirectedGraph> generate_graph_();
    /** assigns "target" symbol to the target vertex and "src<i>" to the i-th child vertex. Also
        appends these symbols to S. */
    void assign_symbols_(SafePtr<CodeSymbols>& S);
    /** given an ImplicitDimension for the computation, adapt it for this recurrence
        relation. Default version does not do anything. */
    virtual SafePtr<ImplicitDimensions> adapt_dims_(const SafePtr<ImplicitDimensions>& dims) const;

    /// does this recurrent relation have a generic equivalent? Default is no.
    virtual bool has_generic(const SafePtr<CompilationParameters>& cparams) const;
    /// return the name of a header file with the declaration of the generic code
    virtual std::string generic_header() const;
    /// return the implementation of this recurrence relation in terms of generic code
    virtual std::string generic_instance(const SafePtr<CodeContext>& context, const SafePtr<CodeSymbols>& args) const;
    
  };

  namespace algebra {
    /// these operators are extremely useful to write compact expressions
    SafePtr<RecurrenceRelation::ExprType> operator+(const SafePtr<DGVertex>& A,
                                                    const SafePtr<DGVertex>& B);
    SafePtr<RecurrenceRelation::ExprType> operator-(const SafePtr<DGVertex>& A,
                                                    const SafePtr<DGVertex>& B);
    SafePtr<RecurrenceRelation::ExprType> operator*(const SafePtr<DGVertex>& A,
                                                    const SafePtr<DGVertex>& B);
    SafePtr<RecurrenceRelation::ExprType> operator/(const SafePtr<DGVertex>& A,
                                                    const SafePtr<DGVertex>& B);
    const SafePtr<RecurrenceRelation::ExprType>& operator+=(SafePtr<RecurrenceRelation::ExprType>& A,
                                                            const SafePtr<DGVertex>& B);
    const SafePtr<RecurrenceRelation::ExprType>& operator-=(SafePtr<RecurrenceRelation::ExprType>& A,
                                                            const SafePtr<DGVertex>& B);
    const SafePtr<RecurrenceRelation::ExprType>& operator*=(SafePtr<RecurrenceRelation::ExprType>& A,
                                                            const SafePtr<DGVertex>& B);
    const SafePtr<RecurrenceRelation::ExprType>& operator/=(SafePtr<RecurrenceRelation::ExprType>& A,
                                                            const SafePtr<DGVertex>& B);
  };

  /** RRStack implements a stack of RecurrenceRelation's which can only hold
      one instance of a given RR. RecurrenceRelation::label() is used for hashing
    */
  class RRStack : public SingletonStack<RecurrenceRelation,std::string>
  {
    public:
    typedef SingletonStack<RecurrenceRelation,std::string> parent_type;
    typedef parent_type::data_type data_type;
    typedef parent_type::value_type value_type;
    typedef parent_type::iter_type iter_type;
    typedef parent_type::citer_type citer_type;
    typedef parent_type::InstanceID InstanceID;

    /// Obtain the unique Instance of RRStack
    static SafePtr<RRStack>& Instance();
    ~RRStack() {}
    
    /// adds content of rrs to this stack
    void add(const SafePtr<RRStack>& rrs);
    /// removes rr from the stack
    void remove(const data_type& rr);
    
    private:
    // private constructor because it's a Singleton
    RRStack() : parent_type(&RecurrenceRelation::label) {}

    static SafePtr<RRStack> rrstack_;
  };
  
};

//#include <vrr_11_twoprep_11.h>
//#include <hrr.h>
//#include <shell_to_ints.h>
//#include <iter.h>

#endif

