


#ifndef _libint2_src_bin_libint_dg_h_
#define _libint2_src_bin_libint_dg_h_

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <rr.h>

using namespace std;


namespace libint2 {

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

  class DirectedGraph {

    vector<DGVertex*> stack_;

    static const unsigned int default_size_ = 100;
    unsigned int first_free_;

    // adds a vertex to the graph
    void add_vertex(DGVertex*);
    /** This function is used to implement (recursive) append_target().
        vertex is appended to the graph and then RR is applied to is.
     */
    template <class I, class RR> void recurse(I* vertex);
    /** This function is used to implement (recursive) apply_to_all().
        RR is applied to vertex and all its children.
     */
    template <class RR> void recurse(DGVertex* vertex);

    // Which vertex is the first to compute
    DGVertex* first_to_compute_;
    // prepare_to_traverse must be called before actual traversal
    void prepare_to_traverse();
    // traverse_from(arc) build recurively the traversal order
    void traverse_from(DGArc*);
    // schedule_computation(vertex) puts vertex first in the computation order
    void schedule_computation(DGVertex*);

  public:
    /** This constructor doesn't do much. Actual initialization of the graph
        must be done using append_target */
    DirectedGraph();
    ~DirectedGraph();

    /** non-template append_target appends the vertex to the graph as a target
    */
    void append_target(DGVertex*);

    /** append_target appends I to the graph as a target vertex and applies
        RR to it. append_target can be called multiple times on the same
        graph if more than one target vertex is needed.

        I must derive from DGVertex. RR must derive from RecurrenceRelation.
        RR has a constructor which takes const I& as the only argument.
        RR must have a public member const I* child(unsigned int) .

        NOTE TO SELF : need to implement these restrictions using
        standard Bjarne Stroustrup's approach.

    */
    template <class I, class RR> void append_target(I*);

    /** apply_to_all applies RR to all vertices already on the graph.

        RR must derive from RecurrenceRelation. RR must define TargetType
        as a typedef.
        RR must have a public member const DGVertex* child(unsigned int) .

        NOTE TO SELF : need to implement these restrictions using
        standard Bjarne Stroustrup's approach.

    */
    template <class RR> void apply_to_all();

    /** after all append_target's have been called, traverse()
        construct a heuristic order of traversal for the graph.

    */
    DGVertex* traverse();

    /// Prints out call sequence
    void debug_print_traversal(ostream& os) const;

    /// Resets the graph and all vertices
    void reset();
  };

  /// Apply RR to target
  template <class I, class RR>
    void
    DirectedGraph::append_target(I* target)
    {
      target->make_a_target();
      recurse<I,RR>(target);
    };

  /// Apply RR to target
  template <class I, class RR>
    void
    DirectedGraph::recurse(I* vertex)
    {
      try {
        add_vertex(vertex);
      }
      catch (VertexAlreadyOnStack) {
        return;
      }
      
      RR* rr0 = new RR(vertex);
      const int num_children = rr0->num_children();
      
      for(int c=0; c<num_children; c++) {
        
        DGVertex* child = rr0->child(c);
        DGArc* arc = new DGArcRel<RR>(vertex,child,rr0);
        vertex->add_exit_arc(arc);
        
        recurse<I,RR>(static_cast<I*>(child));
        
      }
    };

  /// Apply RR recursively starting with vertex
  template <class RR>
    void
    DirectedGraph::recurse(DGVertex* vertex)
    {
      try {
        add_vertex(vertex);
      }
      catch (VertexAlreadyOnStack) {
        return;
      }
      
      typedef typename RR::TargetType TT;
      TT* tptr = dynamic_cast< TT* >(vertex);
      if (tptr == 0)
        return;
      
      RR* rr0 = new RR(tptr);
      const int num_children = rr0->num_children();
      
      for(int c=0; c<num_children; c++) {
        
        DGVertex* child = rr0->child(c);
        DGArc* arc = new DGArcRel<RR>(vertex,child,rr0);
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
      for(int v=0; v<first_free_; v++) {
        if (stack_[v]->num_exit_arcs() != 0)
          continue;
        TT* tptr = dynamic_cast< TT* >(stack_[v]);
        if (tptr == 0)
          continue;
      
        RR* rr0 = new RR(tptr);
        const int num_children = rr0->num_children();
      
        for(int c=0; c<num_children; c++) {
        
          DGVertex* child = rr0->child(c);
          DGArc* arc = new DGArcRel<RR>(tptr,child,rr0);
          tptr->add_exit_arc(arc);
        
          recurse<RR>(child);
        
        }
      }
    }

};


#endif
