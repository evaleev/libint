/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint compiler.
 *
 *  Libint compiler is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint compiler is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint compiler.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_bin_libint_dg_h_
#define _libint2_src_bin_libint_dg_h_

#include <dgvertex.h>
#include <exception.h>
#include <global_macros.h>
#include <key.h>
#include <smart_ptr.h>

#include <algorithm>
#include <cassert>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace libint2 {

//  class DGVertex;
class DGArc;
template <class T>
class DGArcRel;
template <class T>
class AlgebraicOperator;
class Strategy;
class Tactic;
class CodeContext;
class MemoryManager;
class ImplicitDimensions;
class CodeSymbols;
class GraphRegistry;
class InternalGraphRegistry;

/** DirectedGraph is an implementation of a directed graph
    composed of vertices represented by DGVertex objects. Most important
   operations will assume that this is a DAG, i.e. there are no directed cycles.

    \note The objects are allocated on free store and the graph is implemented
   as an object of type 'vertices'.
 */

class DirectedGraph : public std::enable_shared_from_this<DirectedGraph> {
 public:
  typedef DGVertex vertex;
  typedef DGArc arc;
  typedef std::shared_ptr<DGVertex> ver_ptr;
  typedef std::shared_ptr<DGArc> arc_ptr;
  typedef DGVertexKey key_type;
  typedef std::multimap<key_type, ver_ptr> VPtrAssociativeContainer;
  typedef std::list<ver_ptr> VPtrSequenceContainer;

  typedef VPtrSequenceContainer targets;
#if USE_ASSOCCONTAINER_BASED_DIRECTEDGRAPH
  typedef VPtrAssociativeContainer vertices;
#else
  typedef VPtrSequenceContainer vertices;
#endif
  typedef targets::iterator target_iter;
  typedef targets::const_iterator target_citer;
  typedef vertices::iterator ver_iter;
  typedef vertices::const_iterator ver_citer;
  // not possible: typedef vertex::Address address;
  typedef int address;
  // not possible: typedef vertex::Size size;
  typedef unsigned int size;
  typedef std::vector<address> addresses;

 private:
  /// converts what is stored in the container to a smart ptr to the vertex
  static inline const ver_ptr& vertex_ptr(
      const VPtrAssociativeContainer::value_type& v) {
    return v.second;
  }
  static inline const ver_ptr& vertex_ptr(
      const VPtrSequenceContainer::value_type& v) {
    return v;
  }
  /// converts what is stored in the container to a smart ptr to the vertex
  static inline ver_ptr& vertex_ptr(VPtrAssociativeContainer::value_type& v) {
    return v.second;
  }
  static inline ver_ptr& vertex_ptr(VPtrSequenceContainer::value_type& v) {
    return v;
  }

 public:
  /** Creates an empty DAG. Actual initialization of the graph
      must be done using append_target */
  DirectedGraph();
  ~DirectedGraph();

  /// Returns the number of vertices
  unsigned int num_vertices() const { return stack_.size(); }
#if 0
    /// Returns all vertices
    const vertices& all_vertices() const { return stack_; }
    /// Returns all targets
    const targets& all_targets() const { return targets_; }
#endif
  /// Find vertex v or it's equivalent on the graph. Return null pointer if not
  /// found. Computational cost for a graph based on a nonassociative container
  /// may be high
  const std::shared_ptr<DGVertex>& find(
      const std::shared_ptr<DGVertex>& v) const {
    return vertex_is_on(v);
  }

  /** appends v to the graph. If v's copy is already on the graph, return the
     pointer to the copy. Else return pointer to *v.
   */
  std::shared_ptr<DGVertex> append_vertex(const std::shared_ptr<DGVertex>& v);

  /** non-template append_target appends the vertex to the graph as a target
   */
  void append_target(const std::shared_ptr<DGVertex>&);

  /** append_target appends I to the graph as a target vertex and applies
      RR to it. append_target can be called multiple times on the same
      graph if more than one target vertex is needed.

      I must derive from DGVertex. RR must derive from RecurrenceRelation.
      RR has a constructor which takes const I& as the only argument.
      RR must have a public member const I* child(unsigned int) .

      NOTE TO SELF : need to implement these restrictions using
      standard Bjarne Stroustrup's approach.

   */
  template <class I, class RR>
  void append_target(const std::shared_ptr<I>& target) {
    target->make_a_target();
    recurse<I, RR>(target);
  }

  /** apply_to_all applies RR to all vertices already on the graph.

      RR must derive from RecurrenceRelation. RR must define TargetType
      as a typedef.
      RR must have a public member const DGVertex* child(unsigned int) .

      NOTE TO SELF : need to implement these restrictions using
      standard Bjarne Stroustrup's approach.

   */
  template <class RR>
  void apply_to_all() {
    typedef typename RR::TargetType TT;
    typedef vertices::iterator iter;
    for (iter v = stack_.begin(); v != stack_.end(); ++v) {
      ver_ptr& vptr = vertex_ptr(*v);
      if ((vptr)->num_exit_arcs() != 0) continue;
      std::shared_ptr<TT> tptr = std::dynamic_pointer_cast<TT, DGVertex>(v);
      if (tptr == 0) continue;

      std::shared_ptr<RR> rr0(new RR(tptr));
      const int num_children = rr0->num_children();

      for (int c = 0; c < num_children; c++) {
        std::shared_ptr<DGVertex> child = rr0->child(c);
        std::shared_ptr<DGArc> arc(new DGArcRel<RR>(tptr, child, rr0));
        tptr->add_exit_arc(arc);

        recurse<RR>(child);
      }
    }
  }

  /** after all append_target's have been called, apply(strategy,tactic)
    constructs a graph. strategy specifies how to apply recurrence relations.
    The goal of strategies is to connect the target vertices to simpler,
    precomputable vertices. There usually are many ways to reduce a vertex.
    Tactic specifies which of these possibilities to choose.
   */
  void apply(const std::shared_ptr<Strategy>& strategy,
             const std::shared_ptr<Tactic>& tactic);

  typedef void (DGVertex::*DGVertexMethodPtr)();
  /** apply_at<method>(vertex) calls method() on vertex and all of its
   * descendants
   */
  template <DGVertexMethodPtr method>
  void apply_at(const std::shared_ptr<DGVertex>& vertex) const {
    ((vertex.get())->*method)();
    typedef DGVertex::ArcSetType::const_iterator aciter;
    const aciter abegin = vertex->first_exit_arc();
    const aciter aend = vertex->plast_exit_arc();
    for (aciter a = abegin; a != aend; ++a) apply_at<method>((*a)->dest());
  }

  /** calls Method(v) for each v, iterating in forward direction */
  template <class Method>
  void foreach (Method& m) {
    typedef vertices::iterator iter;
    for (iter v = stack_.begin(); v != stack_.end(); ++v) {
      ver_ptr& vptr = vertex_ptr(*v);
      m(vptr);
    }
  }

  /** calls Method(v) for each v, iterating in forward direction */
  template <class Method>
  void foreach (Method& m) const {
    typedef vertices::const_iterator citer;
    for (citer v = stack_.begin(); v != stack_.end(); ++v) {
      const ver_ptr& vptr = vertex_ptr(*v);
      m(vptr);
    }
  }
  /** calls Method(v) for each v, iterating in reverse direction */
  template <class Method>
  void rforeach(Method& m) {
    typedef vertices::reverse_iterator riter;
    for (riter v = stack_.rbegin(); v != stack_.rend(); ++v) {
      ver_ptr& vptr = vertex_ptr(*v);
      m(vptr);
    }
  }

  /** calls Method(v) for each v, iterating in reverse direction */
  template <class Method>
  void rforeach(Method& m) const {
    typedef vertices::const_reverse_iterator criter;
    for (criter v = stack_.rbegin(); v != stack_.rend(); ++v) {
      const ver_ptr& vptr = vertex_ptr(*v);
      m(vptr);
    }
  }

  /** after Strategy has been applied, simple recurrence relations need to be
      optimized away. optimize_rr_out() will replace all simple recurrence
     relations with code representing them.
   */
  void optimize_rr_out(const std::shared_ptr<CodeContext>& context);

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
  void generate_code(const std::shared_ptr<CodeContext>& context,
                     const std::shared_ptr<MemoryManager>& memman,
                     const std::shared_ptr<ImplicitDimensions>& dims,
                     const std::shared_ptr<CodeSymbols>& args,
                     const std::string& label, std::ostream& decl,
                     std::ostream& code);

  /** Resets the graph and all vertices. The stack of unresolved recurrence
      relations is preserved.
   */
  void reset();

  /** num_children_on(rr) returns the number of children of rr which
      are already on this graph.
   */
  template <class RR>
  unsigned int num_children_on(const std::shared_ptr<RR>& rr) const {
    unsigned int nchildren = rr->num_children();
    unsigned int nchildren_on_stack = 0;
    for (unsigned int c = 0; c < nchildren; c++) {
      if (!vertex_is_on(rr->rr_child(c)))
        continue;
      else
        nchildren_on_stack++;
    }

    return nchildren_on_stack;
  }

  /// Returns the registry
  std::shared_ptr<GraphRegistry>& registry() { return registry_; }
  const std::shared_ptr<GraphRegistry>& registry() const { return registry_; }

  /// return the graph label
  const std::string& label() const { return label_; }
  /// sets label to \c new_label
  void set_label(const std::string& new_label) { label_ = new_label; }

  /// return true if there are vertices with 0 children but not pre-computed
  bool missing_prerequisites() const;

 private:
  /// contains vertices
  vertices stack_;
  /// refers to targets, cannot be an associative container -- order of
  /// iteration over targets is important
  targets targets_;
  /// addresses of blocks which accumulate targets
  addresses target_accums_;

  // graph label, used for annotating internal work, e.g. graphviz plots
  std::string label_;

  typedef std::map<std::string, bool> FuncNameContainer;
  /** Maintains the list of names of functions calls to which have been
     generated so far. It is used to generate include statements.
   */
  FuncNameContainer func_names_;

#if !USE_ASSOCCONTAINER_BASED_DIRECTEDGRAPH
  static const unsigned int default_size_ = 100;
#endif

  // maintains data about the graph which does not belong IN the graph
  std::shared_ptr<GraphRegistry> registry_;
  // maintains private data about the graph which does not belong IN the graph
  std::shared_ptr<InternalGraphRegistry> iregistry_;

  /// Access the internal registry
  std::shared_ptr<InternalGraphRegistry>& iregistry() { return iregistry_; }
  const std::shared_ptr<InternalGraphRegistry>& iregistry() const {
    return iregistry_;
  }

  /** adds a vertex to the graph. If the vertex already found on the graph
      then the vertex is not added and the function returns false */
  std::shared_ptr<DGVertex> add_vertex(const std::shared_ptr<DGVertex>&);
  /** same as add_vertex(), only assumes that there's no equivalent vertex on
   * the graph (see vertex_is_on) */
  void add_new_vertex(const std::shared_ptr<DGVertex>&);
  /// returns true if vertex if already on graph
  const std::shared_ptr<DGVertex>& vertex_is_on(
      const std::shared_ptr<DGVertex>& vertex) const;
  /// removes vertex from the graph. may throw CannotPerformOperation
  void del_vertex(vertices::iterator&);
  /** This function is used to implement (recursive) append_target().
      vertex is appended to the graph and then RR is applied to is.
   */
  template <class I, class RR>
  void recurse(const std::shared_ptr<I>& vertex) {
    std::shared_ptr<DGVertex> dgvertex = add_vertex(vertex);
    if (dgvertex != vertex) return;

    std::shared_ptr<RR> rr0(new RR(vertex));
    const int num_children = rr0->num_children();

    for (int c = 0; c < num_children; c++) {
      std::shared_ptr<DGVertex> child = rr0->child(c);
      std::shared_ptr<DGArc> arc(new DGArcRel<RR>(vertex, child, rr0));
      vertex->add_exit_arc(arc);

      std::shared_ptr<I> child_cast =
          std::dynamic_pointer_cast<I, DGVertex>(child);
      if (child_cast == 0)
        throw std::runtime_error(
            "DirectedGraph::recurse(const std::shared_ptr<I>& vertex) -- "
            "dynamic cast failed, most probably this is a logic error!");
      recurse<I, RR>(child_cast);
    }
  }

  /** This function is used to implement (recursive) apply_to_all().
      RR is applied to vertex and all its children.
   */
  template <class RR>
  void recurse(const std::shared_ptr<DGVertex>& vertex) {
    std::shared_ptr<DGVertex> dgvertex = add_vertex(vertex);
    if (dgvertex != vertex) return;

    typedef typename RR::TargetType TT;
    std::shared_ptr<TT> tptr = std::dynamic_pointer_cast<TT, DGVertex>(vertex);
    if (tptr == 0) return;

    std::shared_ptr<RR> rr0(new RR(tptr));
    const int num_children = rr0->num_children();

    for (int c = 0; c < num_children; c++) {
      std::shared_ptr<DGVertex> child = rr0->child(c);
      std::shared_ptr<DGArc> arc(new DGArcRel<RR>(vertex, child, rr0));
      vertex->add_exit_arc(arc);

      recurse<RR>(child);
    }
  }

  /** This function is used to implement (recursive) apply().
    strategy and tactic are applied to vertex and all its children.
   */
  void apply_to(const std::shared_ptr<DGVertex>& vertex,
                const std::shared_ptr<Strategy>& strategy,
                const std::shared_ptr<Tactic>& tactic);
  /// This function insert expr of type AlgebraicOperator<DGVertex> into the
  /// graph
  std::shared_ptr<DGVertex> insert_expr_at(
      const std::shared_ptr<DGVertex>& where,
      const std::shared_ptr<AlgebraicOperator<DGVertex> >& expr);
  /// This function replaces RecurrenceRelations with concrete arithemtical
  /// expressions
  void replace_rr_with_expr();
  /// This function gets rid of trivial math such as multiplication/division
  /// by 1.0, etc.
  void remove_trivial_arithmetics();
  /** This function gets rid of nodes which are connected
  to their equivalents (such as (ss|ss) shell quartet can only be connected to
  (ss|ss) integral)
   */
  void handle_trivial_nodes(const std::shared_ptr<CodeContext>& context);
  /// This functions removes vertices not connected to other vertices
  void remove_disconnected_vertices();
  /** Finds (binary) subtrees. The subtrees correspond to a single-line code (no
     intermediates are used in other expressions)
   */
  void find_subtrees();
  /** Finds (binary) subtrees starting (recursively) at v.
   */
  void find_subtrees_from(const std::shared_ptr<DGVertex>& v);
  /** If v1 and v2 are connected by DGArcDirect and all entry arcs to v1 are of
     the DGArcDirect type as well, this function will reattach all arcs extering
     v1 to v2 and remove v1 from the graph altogether. May throw
     CannotPerformOperation.
   */
  bool remove_vertex_at(const std::shared_ptr<DGVertex>& v1,
                        const std::shared_ptr<DGVertex>& v2);

  // Which vertex is the first to compute
  std::shared_ptr<DGVertex> first_to_compute_;
  // prepare_to_traverse must be called before actual traversal
  void prepare_to_traverse();
  // traverse_from(arc) build recurively the traversal order
  void traverse_from(const std::shared_ptr<DGArc>&);
  // schedule_computation(vertex) puts vertex first in the computation order
  void schedule_computation(const std::shared_ptr<DGVertex>&);

  // Compute addresses on stack assuming that quantities larger than
  // min_size_to_alloc to be allocated on stack
  void allocate_mem(const std::shared_ptr<MemoryManager>& memman,
                    const std::shared_ptr<ImplicitDimensions>& dims,
                    unsigned int min_size_to_alloc = 1);
  // Assign symbols to the vertices
  void assign_symbols(const std::shared_ptr<CodeContext>& context,
                      const std::shared_ptr<ImplicitDimensions>& dims);
  // If v is an AlgebraicOperator, assign (recursively) symbol to the operator.
  // All other must have been already assigned
  void assign_oper_symbol(const std::shared_ptr<CodeContext>& context,
                          std::shared_ptr<DGVertex>& v);
  // Print the code using symbols generated with assign_symbols()
  void print_def(const std::shared_ptr<CodeContext>& context, std::ostream& os,
                 const std::shared_ptr<ImplicitDimensions>& dims,
                 const std::shared_ptr<CodeSymbols>& args);

  /** Returns true if cannot enclose the code in a vector loop
      Possible reason: the traversal path contains a RecurrenceRelation that
     generates a function call (most do except IntegralSet_to_Integrals; \sa
     RecurrentRelation::is_simple() )
      */
  bool cannot_enclose_in_outer_vloop() const;
};

//
// Nonmember utilities
//

/// converts what is stored in the container to a smart ptr to the vertex
inline const DirectedGraph::ver_ptr& vertex_ptr(
    const DirectedGraph::VPtrAssociativeContainer::value_type& v) {
  return v.second;
}
inline const DirectedGraph::ver_ptr& vertex_ptr(
    const DirectedGraph::VPtrSequenceContainer::value_type& v) {
  return v;
}
/// converts what is stored in the container to a smart ptr to the vertex
inline DirectedGraph::ver_ptr& vertex_ptr(
    DirectedGraph::VPtrAssociativeContainer::value_type& v) {
  return v.second;
}
inline DirectedGraph::ver_ptr& vertex_ptr(
    DirectedGraph::VPtrSequenceContainer::value_type& v) {
  return v;
}

#if USE_ASSOCCONTAINER_BASED_DIRECTEDGRAPH
inline DirectedGraph::key_type key(const DGVertex& v);
#endif

//
// Nonmember predicates
//

/// return true if there are non-unrolled targets
bool nonunrolled_targets(const DirectedGraph::targets& targets);

/// extracts external symbols and RRs from the graph
void extract_symbols(const std::shared_ptr<DirectedGraph>& dg);

// use these functors with DirectedGraph::foreach
struct PrerequisitesExtractor {
  std::deque<std::shared_ptr<DGVertex> > vertices;
  void operator()(const std::shared_ptr<DGVertex>& v);
};
struct VertexPrinter {
  VertexPrinter(std::ostream& ostr) : os(ostr) {}
  std::ostream& os;
  void operator()(const std::shared_ptr<DGVertex>& v);
};

};  // namespace libint2

#endif
