
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <assert.h>
#include <exception.h>
#include <smart_ptr.h>

#ifndef _libint2_src_bin_libint_dg_h_
#define _libint2_src_bin_libint_dg_h_

using namespace std;


namespace libint2 {

  class DGVertex;
  class DGArc;
  template <class T> class DGArcRel;
  template <class T> class AlgebraicOperator;
  class Strategy;
  class Tactic;
  class CodeContext;
  class MemoryManager;
  class ImplicitDimensions;
  class CodeSymbols;
  class GraphRegistry;

  /** DirectedGraph is an implementation of a directed graph
      composed of vertices represented by DGVertex objects. The objects
      are allocated on free store and the graph is implemented as
      an object of type container.
  */

  class DirectedGraph : public EnableSafePtrFromThis<DirectedGraph> {
  public:
    typedef DGVertex vertex;
    typedef DGArc arc;
    typedef SafePtr<DGVertex> ver_ptr;
    typedef SafePtr<DGArc> arc_ptr;
    typedef vector<ver_ptr> container;
    typedef container::iterator iter;
    typedef container::const_iterator citer;
    
    /** This constructor doesn't do much. Actual initialization of the graph
        must be done using append_target */
    DirectedGraph();
    ~DirectedGraph();

    /// Returns the number of vertices
    const unsigned int num_vertices() const { return first_free_; }
    
    /** appends v to the graph
    */
    void append_vertex(const SafePtr<DGVertex>& v) throw(VertexAlreadyOnStack);

    /** non-template append_target appends the vertex to the graph as a target
    */
    void append_target(const SafePtr<DGVertex>&);

    /** append_target appends I to the graph as a target vertex and applies
        RR to it. append_target can be called multiple times on the same
        graph if more than one target vertex is needed.

        I must derive from DGVertex. RR must derive from RecurrenceRelation.
        RR has a constructor which takes const I& as the only argument.
        RR must have a public member const I* child(unsigned int) .

        NOTE TO SELF : need to implement these restrictions using
        standard Bjarne Stroustrup's approach.

    */
    template <class I, class RR> void append_target(const SafePtr<I>&);

    /** apply_to_all applies RR to all vertices already on the graph.

        RR must derive from RecurrenceRelation. RR must define TargetType
        as a typedef.
        RR must have a public member const DGVertex* child(unsigned int) .

        NOTE TO SELF : need to implement these restrictions using
        standard Bjarne Stroustrup's approach.

    */
    template <class RR> void apply_to_all();
    
    /** after all append_target's have been called, apply(strategy,tactic)
      constructs a graph. strategy specifies how to apply recurrence relations.
      The goal of strategies is to connect the target vertices to simpler, precomputable vertices.
      There usually are many ways to reduce a vertex.
      Tactic specifies which of these possibilities to choose.
      */
    void apply(const SafePtr<Strategy>& strategy,
               const SafePtr<Tactic>& tactic);
    
    typedef void (DGVertex::* DGVertexMethodPtr)();
    /** apply_at<method>(vertex) calls method() on vertex and all of its descendants
      */
    template <DGVertexMethodPtr method>
      void apply_at(const SafePtr<DGVertex>& vertex) const;

    /** after Strategy has been applied, simple recurrence relations need to be
        optimized away. optimize_rr_out() will replace all simple recurrence relations
        with code representing them.
    */
    void optimize_rr_out();

    /** after all apply's have been called, traverse()
        construct a heuristic order of traversal for the graph.
    */
    void traverse();
    
    /** update func_names_
    */
    void update_func_names();

    /// Prints out call sequence
    void debug_print_traversal(ostream& os) const;
    
    /**
    Prints out the graph in format understood by program "dot"
    of package "graphviz". If symbols is true then label vertices
    using their symbols rather than (descriptive) labels.
    */
    void print_to_dot(bool symbols, std::ostream& os = std::cout) const;
    
    /**
       Generates code for the current computation using context.
       dims specifies the implicit dimensions,
       args specifies the code symbols for the arguments to the function,
       label specifies the tag for the computation,
       decl specifies the stream to receive declaration code,
       code specifies the stream to receive the definition code
    */
    void generate_code(const SafePtr<CodeContext>& context, const SafePtr<MemoryManager>& memman,
                       const SafePtr<ImplicitDimensions>& dims, const SafePtr<CodeSymbols>& args,
                       const std::string& label,
                       std::ostream& decl, std::ostream& code);
    
    /** Resets the graph and all vertices. The stack of unresolved recurrence
        relations is preserved.
      */
    void reset();

    /** num_children_on(rr) returns the number of children of rr which
        are already on this graph.
    */
    template <class RR>
      unsigned int
      num_children_on(const SafePtr<RR>& rr) const;
    
    /// Returns the registry
    SafePtr<GraphRegistry>& registry() { return registry_; }
    const SafePtr<GraphRegistry>& registry() const { return registry_; }
    
  private:

    /// contains vertices
    container stack_;

    typedef std::map<std::string,bool> FuncNameContainer;
    /** Maintains the list of names of functions calls to which have been generated so far.
        It is used to generate include statements.
    */
    FuncNameContainer func_names_;
    
    static const unsigned int default_size_ = 100;
    unsigned int first_free_;
    
    // maintains data about the graph which does not belong IN the graph
    SafePtr<GraphRegistry> registry_;
    
    /** adds a vertex to the graph. If the vertex already found on the graph
        then the vertex is not added and the function returns false */
    void add_vertex(const SafePtr<DGVertex>&) throw(VertexAlreadyOnStack);
    /** same as add_vertex(), only assumes that there's no equivalent vertex on the graph (see vertex_is_on) */
    void add_new_vertex(const SafePtr<DGVertex>&);
    /// returns true if vertex if already on graph
    void vertex_is_on(const SafePtr<DGVertex>& vertex) const throw(VertexAlreadyOnStack);
    /// removes vertex from the graph
    void del_vertex(const SafePtr<DGVertex>&) throw(CannotPerformOperation);
    /** This function is used to implement (recursive) append_target().
        vertex is appended to the graph and then RR is applied to is.
     */
    template <class I, class RR> void recurse(const SafePtr<I>& vertex);
    /** This function is used to implement (recursive) apply_to_all().
        RR is applied to vertex and all its children.
     */
    template <class RR> void recurse(const SafePtr<DGVertex>& vertex);
    /** This function is used to implement (recursive) apply().
      strategy and tactic are applied to vertex and all its children.
    */
    void apply_to(const SafePtr<DGVertex>& vertex,
                  const SafePtr<Strategy>& strategy,
                  const SafePtr<Tactic>& tactic);
    /// This function insert expr of type AlgebraicOperator<DGVertex> into the graph
    void insert_expr_at(const SafePtr<DGVertex>& where, const SafePtr< AlgebraicOperator<DGVertex> >& expr) throw(VertexAlreadyOnStack);
    /// This function replaces RecurrenceRelations with concrete arithemtical expressions
    void replace_rr_with_expr();
    /// This function gets rid of trivial math such as multiplication/division by 1.0, etc.
    void remove_trivial_arithmetics();
    /** This function gets rid of nodes which are connected
    to their equivalents (such as (ss|ss) shell quartet can only be connected to (ss|ss) integral)
    */
    void handle_trivial_nodes();
    /// This functions removes vertices not connected to other vertices
    void remove_disconnected_vertices();
    /** If v1 and v2 are connected by DGArcDirect and all entry arcs to v1 are of the DGArcDirect type as well,
        this function will reattach all arcs extering v1 to v2 and remove v1 from the graph alltogether. */
    void remove_vertex_at(const SafePtr<DGVertex>& v1, const SafePtr<DGVertex>& v2) throw(CannotPerformOperation);
    
    // Which vertex is the first to compute
    SafePtr<DGVertex> first_to_compute_;
    // prepare_to_traverse must be called before actual traversal
    void prepare_to_traverse();
    // traverse_from(arc) build recurively the traversal order
    void traverse_from(const SafePtr<DGArc>&);
    // schedule_computation(vertex) puts vertex first in the computation order
    void schedule_computation(const SafePtr<DGVertex>&);

    // Compute addresses on stack assuming that quantities larger than min_size_to_alloc to be allocated on stack
    void allocate_mem(const SafePtr<MemoryManager>& memman,
                      const SafePtr<ImplicitDimensions>& dims,
                      unsigned int min_size_to_alloc = 1);
    // Assign symbols to the vertices
    void assign_symbols(const SafePtr<CodeContext>& context, const SafePtr<ImplicitDimensions>& dims);
    // Print the code using symbols generated with assign_symbols()
    void print_def(const SafePtr<CodeContext>& context, std::ostream& os,
                   const SafePtr<ImplicitDimensions>& dims);
    
    // Returns true if the traversal path contains a nontrivial RecurrenceRelation (i.e. not of IntegralSet_to_Integrals variety)
    bool contains_nontrivial_rr() const;

  };

};


#endif
