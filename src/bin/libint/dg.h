
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <smart_ptr.h>
#include <rr.h>

#ifndef _libint2_src_bin_libint_dg_h_
#define _libint2_src_bin_libint_dg_h_

using namespace std;


namespace libint2 {

  class Strategy;

  class CannotAddArc : public std::logic_error {
    
    public:
    CannotAddArc(const std::string& a) :
      logic_error(a) {};
    
  };

  class VertexAlreadyOnStack : public std::logic_error {
    
    public:
    VertexAlreadyOnStack(const std::string& a) :
      logic_error(a) {};
    
  };

  /** DirectedGraph is an implementation of a directed graph
      composed of vertices represented by DGVertex objects. The objects
      are allocated on free store and the graph is implemented as
      vector<DGVertex*>.
  */

  class DirectedGraph : public EnableSafePtrFromThis<DirectedGraph> {

    vector< SafePtr<DGVertex> > stack_;

    static const unsigned int default_size_ = 100;
    unsigned int first_free_;

    /** adds a vertex to the graph. If the vertex already found on the graph
        then the vertex is not added and the routine returns false */
    bool add_vertex(const SafePtr<DGVertex>&);
    // returns true if vertex if already on graph
    bool vertex_is_on(const SafePtr<DGVertex>& vertex) const;
    /** This function is used to implement (recursive) append_target().
        vertex is appended to the graph and then RR is applied to is.
     */
    template <class I, class RR> void recurse(const SafePtr<I>& vertex);
    /** This function is used to implement (recursive) apply_to_all().
        RR is applied to vertex and all its children.
     */
    template <class RR> void recurse(const SafePtr<DGVertex>& vertex);
    /** This function is used to implement (recursive) apply().
      strategy is applied to vertex and all its children.
    */
    void apply_to(const SafePtr<DGVertex>& vertex, const SafePtr<Strategy>& strategy);
    /// This function insert expr of type AlgebraicOperator<DGVertex> into the graph
    void insert_expr_at(const SafePtr<DGVertex>& where, const SafePtr< AlgebraicOperator<DGVertex> >& expr);
    
    // Which vertex is the first to compute
    SafePtr<DGVertex> first_to_compute_;
    // prepare_to_traverse must be called before actual traversal
    void prepare_to_traverse();
    // traverse_from(arc) build recurively the traversal order
    void traverse_from(const SafePtr<DGArc>&);
    // schedule_computation(vertex) puts vertex first in the computation order
    void schedule_computation(const SafePtr<DGVertex>&);

  public:
    /** This constructor doesn't do much. Actual initialization of the graph
        must be done using append_target */
    DirectedGraph();
    ~DirectedGraph();

    /// Returns the number of vertices
    const unsigned int num_vertices() const { return first_free_; }

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
    
    /** after all append_target's have been called, apply(const SafePtr<Strategy>&)
      constructs a graph using Strategy. Strategy specifies how to apply recurrence relations.
      The goal of strategies is to connect the target vertices to simpler, precomputable vertices.
      */
    void apply(const SafePtr<Strategy>&);

    /** after Strategy has been applied, simple recurrence relations need to be
        optimized away. optimize_rr_out() will replace all simple recurrence relations
        with code representing them.
    */
    void optimize_rr_out();

    /** after all apply's have been called, traverse()
        construct a heuristic order of traversal for the graph.
    */
    void traverse();

    /// Prints out call sequence
    void debug_print_traversal(ostream& os) const;
    
    /**
    Prints out the graph in format understood by program "dot"
    of package "graphviz"
    */
    void print_to_dot(ostream& os) const;

    /// Resets the graph and all vertices
    void reset();

    /** num_children_on(rr) returns the number of children of rr which
        are already on this graph.
    */
    template <class RR>
      unsigned int
      num_children_on(const SafePtr<RR>& rr) const;
  };

  /// Apply RR to target
  template <class I, class RR>
    void
    DirectedGraph::append_target(const SafePtr<I>& target)
    {
      target->make_a_target();
      recurse<I,RR>(target);
    };

  /// Apply RR to target
  template <class I, class RR>
    void
    DirectedGraph::recurse(const SafePtr<I>& vertex)
    {
      try {
        add_vertex(vertex);
      }
      catch (VertexAlreadyOnStack) {
        return;
      }
      
      SafePtr<RR> rr0(new RR(vertex));
      const int num_children = rr0->num_children();
      
      for(int c=0; c<num_children; c++) {
        
        SafePtr<DGVertex> child = rr0->child(c);
        SafePtr<DGArc> arc(new DGArcRel<RR>(vertex,child,rr0));
        vertex->add_exit_arc(arc);
        
	SafePtr<I> child_cast = dynamic_pointer_cast<I,DGVertex>(child);
	if (child_cast == 0)
	  throw std::runtime_error("DirectedGraph::recurse(const SafePtr<I>& vertex) -- dynamic cast failed, most probably this is a logic error!");
        recurse<I,RR>(child_cast);
        
      }
    };

  /// Apply RR recursively starting with vertex
  template <class RR>
    void
    DirectedGraph::recurse(const SafePtr<DGVertex>& vertex)
    {
      try {
        add_vertex(vertex);
      }
      catch (VertexAlreadyOnStack) {
        return;
      }
      
      typedef typename RR::TargetType TT;
      SafePtr<TT> tptr = dynamic_pointer_cast<TT,DGVertex>(vertex);
      if (tptr == 0)
        return;
      
      SafePtr<RR> rr0(new RR(tptr));
      const int num_children = rr0->num_children();
      
      for(int c=0; c<num_children; c++) {
        
        SafePtr<DGVertex> child = rr0->child(c);
        SafePtr<DGArc> arc(new DGArcRel<RR>(vertex,child,rr0));
        vertex->add_exit_arc(arc);
        
        recurse<RR>(child);
      }
    };

  /// Apply RR to all classes already on the graph
  template <class RR>
    void
    DirectedGraph::apply_to_all()
    {
      typedef typename RR::TargetType TT;
      const int num_vertices_on_graph = first_free_;
      for(int v=0; v<num_vertices_on_graph; v++) {
        if (stack_[v]->num_exit_arcs() != 0)
          continue;
        SafePtr<TT> tptr = dynamic_pointer_cast<TT,DGVertex>(stack_[v]);
        if (tptr == 0)
          continue;
      
        SafePtr<RR> rr0(new RR(tptr));
        const int num_children = rr0->num_children();
      
        for(int c=0; c<num_children; c++) {
        
          SafePtr<DGVertex> child = rr0->child(c);
          SafePtr<DGArc> arc(new DGArcRel<RR>(tptr,child,rr0));
          tptr->add_exit_arc(arc);
        
          recurse<RR>(child);
        
        }
      }
    }

  template <class RR>
    unsigned int
    DirectedGraph::num_children_on(const SafePtr<RR>& rr) const
    {
      unsigned int nchildren = rr->num_children();
      unsigned int nchildren_on_stack = 0;
      for(int c=0; c<nchildren; c++) {
        if ( vertex_is_on(rr->rr_child(c)) )
          nchildren_on_stack++;
      }

      return nchildren_on_stack;
    }

};


#endif
