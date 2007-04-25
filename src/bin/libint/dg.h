
#ifndef _libint2_src_bin_libint_dg_h_
#define _libint2_src_bin_libint_dg_h_

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <stdexcept>
#include <assert.h>
#include <exception.h>
#include <smart_ptr.h>

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
  class InternalGraphRegistry;

  /** DirectedGraph is an implementation of a directed graph
      composed of vertices represented by DGVertex objects. The objects
      are allocated on free store and the graph is implemented as
      an object of type 'vertices'.
  */

  class DirectedGraph : public EnableSafePtrFromThis<DirectedGraph> {
  public:
    typedef DGVertex vertex;
    typedef DGArc arc;
    typedef SafePtr<DGVertex> ver_ptr;
    typedef SafePtr<DGArc> arc_ptr;
#if USE_MULTIMAP_BASED_DIRECTEDGRAPH
    typedef DGVertexKey key_type;
    typedef std::multimap<key_type,ver_ptr> vertices;
#else
    typedef std::list<ver_ptr> vertices;
#endif
    typedef vertices::iterator ver_iter;
    typedef vertices::const_iterator ver_citer;
    //not possible: typedef vertex::Address address;
    typedef int address;
    //not possible: typedef vertex::Size size;
    typedef unsigned int size;
    typedef std::vector<address> addresses;
    
    /** This constructor doesn't do much. Actual initialization of the graph
        must be done using append_target */
    DirectedGraph();
    ~DirectedGraph();

    /// Returns the number of vertices
    const unsigned int num_vertices() const { return stack_.size(); }
#if 0
    /// Returns all vertices
    const vertices& all_vertices() const { return stack_; }
    /// Returns all targets
    const vertices& all_targets() const { return targets_; }
#endif
    /// Find vertex v or it's equivalent on the graph. Return null pointer if not found. Computational cost for a graph based on a nonassociative container may be high
    const SafePtr<DGVertex>& find(const SafePtr<DGVertex>& v) const { return vertex_is_on(v); }
    
    /** appends v to the graph. If v's copy is already on the graph, return the pointer
	to the copy. Else return pointer to *v.
    */
    SafePtr<DGVertex> append_vertex(const SafePtr<DGVertex>& v);

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

    /** calls Method(v) for each v, iterating in forward direction */
    template <class Method>
      void foreach(Method& m);
    /** calls Method(v) for each v, iterating in forward direction */
    template <class Method>
      void foreach(Method& m) const;
    /** calls Method(v) for each v, iterating in reverse direction */
    template <class Method>
      void rforeach(Method& m);
    /** calls Method(v) for each v, iterating in reverse direction */
    template <class Method>
      void rforeach(Method& m) const;

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
    void debug_print_traversal(std::ostream& os) const;
    
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
    vertices stack_;
    /// refers to targets
    vertices targets_;
    /// addresses of blocks which accumulate targets
    addresses target_accums_;

    typedef std::map<std::string,bool> FuncNameContainer;
    /** Maintains the list of names of functions calls to which have been generated so far.
        It is used to generate include statements.
    */
    FuncNameContainer func_names_;

#if USE_MULTIMAP_BASED_DIRECTEDGRAPH
#else
    static const unsigned int default_size_ = 100;
#endif
    
    // maintains data about the graph which does not belong IN the graph
    SafePtr<GraphRegistry> registry_;
    // maintains private data about the graph which does not belong IN the graph
    SafePtr<InternalGraphRegistry> iregistry_;

    /// Access the internal registry
    SafePtr<InternalGraphRegistry>& iregistry() { return iregistry_; }
    const SafePtr<InternalGraphRegistry>& iregistry() const { return iregistry_; }
    
    /** adds a vertex to the graph. If the vertex already found on the graph
        then the vertex is not added and the function returns false */
    SafePtr<DGVertex> add_vertex(const SafePtr<DGVertex>&);
    /** same as add_vertex(), only assumes that there's no equivalent vertex on the graph (see vertex_is_on) */
    void add_new_vertex(const SafePtr<DGVertex>&);
    /// returns true if vertex if already on graph
    const SafePtr<DGVertex>& vertex_is_on(const SafePtr<DGVertex>& vertex) const;
    /// removes vertex from the graph. may throw CannotPerformOperation
    void del_vertex(const SafePtr<DGVertex>&);
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
    SafePtr<DGVertex> insert_expr_at(const SafePtr<DGVertex>& where, const SafePtr< AlgebraicOperator<DGVertex> >& expr);
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
    /** Finds (binary) subtrees. The subtrees correspond to a single-line code (no intermediates
        are used in other expressions)
    */
    void find_subtrees();
    /** Finds (binary) subtrees starting (recursively) at v.
    */
    void find_subtrees_from(const SafePtr<DGVertex>& v);
    /** If v1 and v2 are connected by DGArcDirect and all entry arcs to v1 are of the DGArcDirect type as well,
        this function will reattach all arcs extering v1 to v2 and remove v1 from the graph altogether.
	May throw CannotPerformOperation.
    */
    void remove_vertex_at(const SafePtr<DGVertex>& v1, const SafePtr<DGVertex>& v2);
    
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
    // If v is an AlgebraicOperator, assign (recursively) symbol to the operator. All other must have been already assigned
    void assign_oper_symbol(const SafePtr<CodeContext>& context, SafePtr<DGVertex>& v);
    // Print the code using symbols generated with assign_symbols()
    void print_def(const SafePtr<CodeContext>& context, std::ostream& os,
                   const SafePtr<ImplicitDimensions>& dims);
    
    // Returns true if the traversal path contains a nontrivial RecurrenceRelation (i.e. not of IntegralSet_to_Integrals variety)
    bool contains_nontrivial_rr() const;

  };


  //
  // Nonmember predicates
  //

  /// return true if there are non-unrolled targets
  bool nonunrolled_targets(const DirectedGraph::vertices& targets);

  /// extracts external symbols and RRs from the graph
  void extract_symbols(const SafePtr<DirectedGraph>& dg);

};


#endif
