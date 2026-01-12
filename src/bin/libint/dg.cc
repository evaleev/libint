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

#include <algebra.h>
#include <codeblock.h>
#include <context.h>
#include <default_params.h>
#include <dg.h>
#include <dims.h>
#include <extract.h>
#include <global_macros.h>
#include <graph_registry.h>
#include <intset_to_ints.h>
#include <prefactors.h>
#include <rr.h>
#include <strategy.h>
#include <task.h>
#include <uncontract.h>

#include <cstdio>
#include <fstream>
#include <functional>
#include <utility>

using namespace std;
using namespace libint2;

#define ONLY_CLONE_IF_DIFF 1

//// utils first

namespace {
void push(DirectedGraph::VPtrAssociativeContainer& vertices,
          const DirectedGraph::ver_ptr& v) {
  DirectedGraph::key_type vkey = libint2::key(*v);
  vertices.insert(std::make_pair(vkey, v));
}
void push(DirectedGraph::VPtrSequenceContainer& vertices,
          const DirectedGraph::ver_ptr& v) {
  vertices.push_back(v);
}
}  // namespace

DirectedGraph::DirectedGraph()
    : stack_(),
      targets_(),
      target_accums_(),
      label_("graph"),
      func_names_(),
      registry_(std::shared_ptr<GraphRegistry>(new GraphRegistry)),
      iregistry_(
          std::shared_ptr<InternalGraphRegistry>(new InternalGraphRegistry)),
      first_to_compute_() {
  stack_.clear();
  targets_.clear();
}

DirectedGraph::~DirectedGraph() { reset(); }

void DirectedGraph::append_target(const std::shared_ptr<DGVertex>& target) {
  target->make_a_target();
  append_vertex(target);
  push(targets_, target);
}

std::shared_ptr<DGVertex> DirectedGraph::append_vertex(
    const std::shared_ptr<DGVertex>& vertex) {
  // If this vertex is owned by this graph, return it immediately
  if (vertex->dg() == this) {
#if DEBUG
    std::cout << "append_vertex: vertex " << vertex->label() << " is already on"
              << endl;
#endif
    return vertex;
  }
  auto vcopy_on_graph = add_vertex(vertex);
  // If this is a new vertex -- tell the vertex who its owner is now
  if (vcopy_on_graph == vertex) vertex->dg(this);
  return vcopy_on_graph;
}

std::shared_ptr<DGVertex> DirectedGraph::add_vertex(
    const std::shared_ptr<DGVertex>& vertex) {
  // if vertex is precomputed -- can replicate it without consequences
  // else check if it's on already
  if (!vertex->precomputed()) {
    std::shared_ptr<DGVertex> vcopy_on_graph = vertex_is_on(vertex);
    if (vcopy_on_graph) return vcopy_on_graph;
  }
  add_new_vertex(vertex);
  return vertex;
}

void DirectedGraph::add_new_vertex(const std::shared_ptr<DGVertex>& vertex) {
#if 0
  // Resize if using std::vector
#if !USE_ASSOCCONTAINER_BASED_DIRECTEDGRAPH
  if (num_vertices() == stack_.capacity()) {
    stack_.resize( stack_.capacity() + default_size_ );
#if DEBUG
    cout << "Increased size of DirectedGraph's stack to "
         << stack_.size() << endl;
#endif
  }
#endif
#endif

  char label[80];
  snprintf(label, 80, "vertex%d", num_vertices());
  vertex->set_graph_label(label);
  vertex->dg(this);

  push(stack_, vertex);
#if DEBUG
  cout << "add_new_vertex: added vertex " << vertex->description() << endl;
#endif

  return;
}

const std::shared_ptr<DGVertex>& DirectedGraph::vertex_is_on(
    const std::shared_ptr<DGVertex>& vertex) const {
  if (vertex->dg() == this) return vertex;

  static std::shared_ptr<DGVertex> null_ptr;
#if USE_ASSOCCONTAINER_BASED_DIRECTEDGRAPH
  typedef vertices::const_iterator citer;
  key_type vkey = key(*vertex);
  // find the first elemnt with this key and iterate until vertex is found or
  // key changes
  citer vpos = stack_.find(vkey);
  const citer end = stack_.end();
  if (vpos != end) {
    bool can_find = true;
    while (can_find) {
      if (can_find && (vpos->second)->equiv(vertex)) {
#if DEBUG
        std::cout << "vertex_is_on: " << (vpos->second)->label() << std::endl;
#endif
        return vpos->second;
      }
      ++vpos;
      can_find = (vpos->first == vkey) && (vpos != end);
    }
  }
#else
  typedef vertices::const_reverse_iterator criter;
  typedef vertices::reverse_iterator riter;
  const criter rend = stack_.rend();
  for (criter v = stack_.rbegin(); v != rend; ++v) {
    const ver_ptr& vptr = vertex_ptr(*v);
    if (vertex->equiv(vptr)) {
#if DEBUG
      std::cout << "vertex_is_on: " << (vptr)->label() << std::endl;
#endif
      return vptr;
    }
  }
#endif
#if DEBUG
  std::cout << "vertex_is_on: NOT " << (vertex)->label() << std::endl;
#endif
  return null_ptr;
}

namespace {
struct __reset_dgvertex {
  void operator()(std::shared_ptr<DGVertex>& v) {
    v->reset();
#if DEBUG
    std::cout << "DirectedGraph::reset: will unregister " << v->label()
              << std::endl;
#endif
    // remove this vertex from its SingletonManager
    v->unregister();
  }
};
struct __reset_shared_ptr {
  void operator()(std::shared_ptr<DGVertex>& v) { v.reset(); }
};
}  // namespace

void DirectedGraph::del_vertex(vertices::iterator& v) {
  static __reset_dgvertex rv;
  if (v == stack_.end())
    throw CannotPerformOperation(
        "DirectedGraph::del_vertex() cannot delete vertex");
  ver_ptr& vptr = vertex_ptr(*v);
  // Cannot delete targets. Should I be able to? Probably not
  if (vptr->is_a_target())
    throw CannotPerformOperation(
        "DirectedGraph::del_vertex() cannot delete targets");
  if (vptr->num_exit_arcs() == 0 && vptr->num_entry_arcs() == 0) {
#if DEBUG
    std::cout << "del_vertex: trying to remove " << (vertex_ptr(*v))->label()
              << std::endl;
#endif
    std::shared_ptr<DGVertex> vptr =
        vertex_ptr(*v);  // keep an instance of the pointer to avoid accidental
                         // automatic destruction of the DGVertex object
    stack_.erase(v);
    rv(vptr);
#if DEBUG
    std::cout << "del_vertex: successful vertex removal " << std::endl;
#endif
  } else
    throw CannotPerformOperation(
        "DirectedGraph::del_vertex() cannot delete vertex");
}

namespace {
struct __prepare_to_traverse {
  void operator()(std::shared_ptr<DGVertex>& v) { v->prepare_to_traverse(); }
};

std::vector<DGVertex::ArcSetType::value_type> sort_children_by_nparents(
    DGVertex::ArcSetType::const_iterator begin,
    DGVertex::ArcSetType::const_iterator end) {
  // std::sort works only for containers that support random access
  //    std::vector<DGVertex::ArcSetType::value_type> sorted_children;
  //    for(DGVertex::ArcSetType::const_iterator i=begin; i!=end; ++i)
  //      sorted_children.push_back(*i);
  std::vector<DGVertex::ArcSetType::value_type> sorted_children(begin, end);
  std::sort(sorted_children.begin(), sorted_children.end(),
            [](DGVertex::ArcSetType::value_type a,
               DGVertex::ArcSetType::value_type b) {
              return a->dest()->num_entry_arcs() < b->dest()->num_entry_arcs();
            });
  return sorted_children;
}
}  // namespace

void DirectedGraph::prepare_to_traverse() {
  __prepare_to_traverse __ptt;
  foreach (__ptt)
    ;
}

/**
 * Recursively traverse depth-first, once a node has been tagged by all of its
 * parents schedule its computation
 */
void DirectedGraph::traverse() {
  // Initialization
  prepare_to_traverse();

  // Start at the targets which don't have parents
  typedef vertices::iterator iter;
  for (iter v = stack_.begin(); v != stack_.end(); ++v) {
    const ver_ptr& vptr = vertex_ptr(*v);
    if ((vptr)->is_a_target() && (vptr)->num_entry_arcs() == 0) {
      // First, since this target doesn't have parents we can schedule its
      // computation
      schedule_computation(vptr);

      //
      // traverse the rest of graph starting from each child
      // traversal will start with children with the fewest parents
      // an explanation for this heuristic: I want to compute the children with
      // 1 parent after other children so that I can potentially fold their
      // computation with addition/subtraction to produce FMA and other
      // composite instructions
      //
      {
        // std::sort works only fro containers that support random access
        std::vector<DGVertex::ArcSetType::value_type> sorted_children =
            sort_children_by_nparents(vptr->first_exit_arc(),
                                      vptr->plast_exit_arc());
        for (auto a = sorted_children.begin(); a != sorted_children.end();
             ++a) {
          traverse_from(*a);
        }
      }
    }
  }
}

void DirectedGraph::traverse_from(const std::shared_ptr<DGArc>& arc) {
  std::shared_ptr<DGVertex> orig = arc->orig();
  std::shared_ptr<DGVertex> dest = arc->dest();
#if DEBUG_TRAVERSAL
  std::cout << "traverse_from: orig = " << orig << " dest = " << dest << endl;
#endif
  // no need to compute if precomputed OR has no children
  if (dest->precomputed() || dest->num_exit_arcs() == 0) {
#if DEBUG_TRAVERSAL
    std::cout << "traverse from: dest is precomputed" << std::endl;
#endif
    return;
  }
  // if has been hit by all parents ...
  const unsigned int num_tags = dest->tag();
#if DEBUG_TRAVERSAL
  std::cout << "traverse from: tagged dest, ntags = " << num_tags << std::endl;
#endif
  const unsigned int num_parents = dest->num_entry_arcs();
  if (num_tags == num_parents) {
    // ... check if it is on a subtree ...
    if (std::shared_ptr<DRTree> stree = dest->subtree()) {
      // ... if yes, schedule only if it is a root of a subtree
      if (stree->root() == dest) schedule_computation(dest);
    }
    // else schedule
    else
      schedule_computation(dest);

    {
      // std::sort works only fro containers that support random access
      std::vector<DGVertex::ArcSetType::value_type> sorted_children =
          sort_children_by_nparents(dest->first_exit_arc(),
                                    dest->plast_exit_arc());
      for (auto a = sorted_children.begin(); a != sorted_children.end(); ++a) {
        traverse_from(*a);
      }
    }
  }
}

void DirectedGraph::schedule_computation(
    const std::shared_ptr<DGVertex>& vertex) {
  vertex->set_postcalc(first_to_compute_);
  first_to_compute_ = vertex;
#if DEBUG || DEBUG_TRAVERSAL
  std::cout << "schedule_computation: " << vertex << endl;
  vertex->print(std::cout);
#endif
}

void DirectedGraph::debug_print_traversal(std::ostream& os) const {
  std::shared_ptr<DGVertex> current_vertex = first_to_compute_;

  os << "Debug print of traversal order" << endl;

  do {
    current_vertex->print(os);
    current_vertex = current_vertex->postcalc();
  } while (current_vertex != 0);
}

namespace {

struct __print_vertices_to_dot {
  bool symbols;
  std::ostream& os;
  __print_vertices_to_dot(bool s, std::ostream& o) : symbols(s), os(o) {}
  void operator()(const std::shared_ptr<DGVertex>& v) {
    os << "  " << v->graph_label() << " [ label = \"";
    if (symbols && v->symbol_set())
      os << v->symbol();
    else
      os << v->label();
    os << "\"]" << endl;
  }
};

struct __print_arcs_to_dot {
  std::ostream& os;
  __print_arcs_to_dot(std::ostream& o) : os(o) {}
  void operator()(const std::shared_ptr<DGVertex>& v) {
    typedef DGVertex::ArcSetType::const_iterator aciter;
    const aciter abegin = v->first_exit_arc();
    const aciter aend = v->plast_exit_arc();
    for (aciter a = abegin; a != aend; ++a) {
      std::shared_ptr<DGVertex> dest = (*a)->dest();
      os << "  " << v->graph_label() << " -> " << dest->graph_label() << endl;
    }
  }
};

}  // namespace

void DirectedGraph::print_to_dot(bool symbols, std::ostream& os) const {
  os << "digraph G {" << endl << "  size = \"8,8\"" << endl;

  __print_vertices_to_dot pvtd(symbols, os);
  foreach (pvtd)
    ;

  __print_arcs_to_dot patd(os);
  foreach (patd)
    ;

  // Print traversal order using dotted lines
  std::shared_ptr<DGVertex> current_vertex = first_to_compute_;
  if (current_vertex != 0) {
    do {
      std::shared_ptr<DGVertex> next = current_vertex->postcalc();
      if (current_vertex && next) {
        os << "  " << current_vertex->graph_label() << " -> "
           << next->graph_label() << " [ style = dotted constraint = false ]";
      }
      current_vertex = next;
    } while (current_vertex != 0);
  }

  os << endl << "}" << endl;
}

void DirectedGraph::reset() {
  // Reset each vertex, releasing all arcs
  __reset_dgvertex rv;
  foreach (rv)
    ;
  __reset_shared_ptr rptr;
  foreach (rptr)
    ;

  // if everything went OK then empty out stack_ and targets_
  stack_.clear();
  targets_.clear();
  first_to_compute_.reset();
  func_names_.clear();
}

/// Apply a strategy to all vertices not yet computed (i.e. which do not have
/// exit arcs)
void DirectedGraph::apply(const std::shared_ptr<Strategy>& strategy,
                          const std::shared_ptr<Tactic>& tactic) {
  const std::shared_ptr<DirectedGraph> this_ptr = shared_from_this();
  typedef vertices::iterator iter;
  for (iter v = stack_.begin(); v != stack_.end(); ++v) {
    const ver_ptr& vptr = vertex_ptr(*v);
    if ((vptr)->num_exit_arcs() != 0 || (vptr)->precomputed() ||
        !(vptr)->need_to_compute())
      continue;

    std::shared_ptr<RecurrenceRelation> rr0 =
        strategy->optimal_rr(this_ptr, (vptr), tactic);
    if (rr0 == 0) continue;

    // add children to the graph
    std::shared_ptr<DGVertex> target = rr0->rr_target();
    const int num_children = rr0->num_children();
    for (int c = 0; c < num_children; c++) {
      std::shared_ptr<DGVertex> child = rr0->rr_child(c);
      bool new_vertex = true;
      std::shared_ptr<DGVertex> dgchild = append_vertex(child);
      if (dgchild != child) {
        child = dgchild;
        new_vertex = false;
      }
      std::shared_ptr<DGArc> arc(
          new DGArcRel<RecurrenceRelation>(target, child, rr0));
      target->add_exit_arc(arc);
      if (new_vertex) apply_to(child, strategy, tactic);
    }
  }
}

/// Add vertex to graph and apply a strategy to vertex recursively
void DirectedGraph::apply_to(const std::shared_ptr<DGVertex>& vertex,
                             const std::shared_ptr<Strategy>& strategy,
                             const std::shared_ptr<Tactic>& tactic) {
  bool not_yet_computed = !vertex->precomputed() && vertex->need_to_compute() &&
                          (vertex->num_exit_arcs() == 0);
  if (!not_yet_computed) return;
  std::shared_ptr<RecurrenceRelation> rr0 =
      strategy->optimal_rr(shared_from_this(), vertex, tactic);
  if (rr0 == 0) return;

  std::shared_ptr<DGVertex> target = rr0->rr_target();
  const int num_children = rr0->num_children();
  for (int c = 0; c < num_children; c++) {
    std::shared_ptr<DGVertex> child = rr0->rr_child(c);
    bool new_vertex = true;
    std::shared_ptr<DGVertex> dgchild = append_vertex(child);
    if (dgchild != child) {
      child = dgchild;
      new_vertex = false;
    }
    std::shared_ptr<DGArc> arc(
        new DGArcRel<RecurrenceRelation>(target, child, rr0));
    try {
      target->add_exit_arc(arc);
    } catch (...) {
      std::cout << "failed to use RR " << rr0->label() << std::endl;
      throw;
    }
    if (new_vertex) apply_to(child, strategy, tactic);
  }
}

// Optimize out simple recurrence relations
void DirectedGraph::optimize_rr_out(
    const std::shared_ptr<CodeContext>& context) {
  replace_rr_with_expr();
  remove_trivial_arithmetics();
  handle_trivial_nodes(context);
  remove_disconnected_vertices();
  find_subtrees();
}

// Replace recurrence relations with expressions
void DirectedGraph::replace_rr_with_expr() {
  typedef vertices::iterator iter;
  for (iter v = stack_.begin(); v != stack_.end(); ++v) {
    const ver_ptr& vptr = vertex_ptr(*v);
    if ((vptr)->num_exit_arcs()) {
      std::shared_ptr<DGArc> arc0 = *((vptr)->first_exit_arc());
      std::shared_ptr<DGArcRR> arc0_cast =
          std::dynamic_pointer_cast<DGArcRR, DGArc>(arc0);
      if (arc0_cast == 0) continue;
      std::shared_ptr<RecurrenceRelation> rr = arc0_cast->rr();

      // Optimize if the recurrence relation is simple and the target and
      // children are of the same type
      if (rr->is_simple() && rr->invariant_type()) {
#if DEBUG || DEBUG_RESTRUCTURE
        std::cout << "replace_rr_with_expr: replacing " << rr->label() << endl;
        std::cout << "replace_rr_with_expr:      with "
                  << rr->rr_expr()->description() << endl;
        std::cout << "replace_rr_with_expr: nchildren = " << rr->num_children()
                  << endl;
        for (unsigned int c = 0; c < rr->num_children(); ++c) {
          std::cout << "replace_rr_with_expr: child " << c << " "
                    << rr->rr_child(c)->description() << endl;
        }
#endif

        // Remove arcs connecting this vertex to children
        (vptr)->del_exit_arcs();

        // and instead insert the numerical expression
        std::shared_ptr<RecurrenceRelation::ExprType> rr_expr = rr->rr_expr();
        std::shared_ptr<DGVertex> expr_vertex =
            std::static_pointer_cast<RecurrenceRelation::ExprType, DGVertex>(
                rr_expr);
        expr_vertex = insert_expr_at((vptr), rr_expr);
        std::shared_ptr<DGArc> arc(new DGArcDirect((vptr), expr_vertex));
        (vptr)->add_exit_arc(arc);
      }
    }
  }
}

//
// This function is very tricky at the moment. The operands have to be added
// before the operator such that operands are guaranteed to be on graph before
// the operator. This way ExprType::equiv can simply compare operand pointers
// and the operator type (very cheap operations compared to the fully recursive
// explicit comparison).
//
std::shared_ptr<DGVertex> DirectedGraph::insert_expr_at(
    const std::shared_ptr<DGVertex>& where,
    const std::shared_ptr<RecurrenceRelation::ExprType>& expr) {
#if DEBUG
  cout << "insert_expr_at: " << expr->description() << endl;
#endif

  typedef RecurrenceRelation::ExprType ExprType;
  std::shared_ptr<DGVertex> expr_vertex =
      std::static_pointer_cast<DGVertex, ExprType>(expr);

  // If the expression is already on then return it
  if (expr->dg() == this) return expr_vertex;

  std::shared_ptr<DGVertex> left_oper = expr->left();
  std::shared_ptr<DGVertex> right_oper = expr->right();
  bool new_left = (left_oper->dg() != this);    //
  bool new_right = (right_oper->dg() != this);  //
  bool need_to_clone = false;  // clone expression if (parts of) the expression
                               // was found on the graph

  // See if left operand is also an operator
  const std::shared_ptr<ExprType> left_cast =
      std::dynamic_pointer_cast<ExprType, DGVertex>(left_oper);
  // if yes -- add it to the graph recursively
  if (left_cast) {
    left_oper = insert_expr_at(expr_vertex, left_cast);
  }
  // else add it directly
  else {
    left_oper = append_vertex(left_oper);
  }
#if ONLY_CLONE_IF_DIFF
  if (left_oper != expr->left()) {
#if DEBUG
    std::cout << "insert_expr_at: append(left) != left" << std::endl;
#endif
    need_to_clone = true;
    new_left = false;
  }
#else
#error "ONLY_CLONE_IF_DIFF must be true"
#endif

  // See if right operand is also an operator
  std::shared_ptr<ExprType> right_cast =
      std::dynamic_pointer_cast<ExprType, DGVertex>(right_oper);
  // if yes -- add it to the graph recursively
  if (right_cast) {
    right_oper = insert_expr_at(expr_vertex, right_cast);
  }
  // else add it directly
  else {
    right_oper = append_vertex(right_oper);
  }
#if ONLY_CLONE_IF_DIFF
  if (right_oper != expr->right()) {
#if DEBUG
    std::cout << "insert_expr_at: append(right) != right" << std::endl;
#endif
    need_to_clone = true;
    new_right = false;
  }
#else
#error "ONLY_CLONE_IF_DIFF must be true"
#endif

  if (need_to_clone) {
    std::shared_ptr<ExprType> expr_new(
        new ExprType(expr, left_oper, right_oper));
    expr_vertex = std::static_pointer_cast<DGVertex, ExprType>(expr_new);
#if DEBUG
    int nchildren = expr->num_exit_arcs();
    cout << "insert_expr_at: cloned AlgebraicOperator with "
         << expr->num_exit_arcs() << " children" << endl;
    if (nchildren) {
      cout << "Left:  " << expr->left()->description() << endl;
      cout << "Right: " << expr->right()->description() << endl;
    }
#endif
  }

  std::shared_ptr<DGVertex> dgexpr_vertex = expr_vertex;
  if (new_left || new_right)
    add_new_vertex(expr_vertex);
  else {
    const bool do_cse = registry()->do_cse();
    if (do_cse) {
#if DEBUG
      std::cout << "insert_expr_at: appending vertex "
                << expr_vertex->description() << std::endl;
#endif
      dgexpr_vertex = append_vertex(expr_vertex);
    } else {
      add_new_vertex(expr_vertex);
    }
    if (expr_vertex != dgexpr_vertex) {
      if (new_left || new_right) {
        cout << "Problem detected: AlgebraicOperator is found on the stack but "
                "one of its operands was new"
             << endl;
        cout << expr_vertex->description() << endl;
        cout << dgexpr_vertex->description() << endl;
        throw std::runtime_error(
            "DirectedGraph::insert_expr_at() -- vertex is not new but one of "
            "the operands is");
      }
    }
    expr_vertex = dgexpr_vertex;
  }
  std::shared_ptr<DGArc> left_arc(new DGArcDirect(expr_vertex, left_oper));
  expr_vertex->add_exit_arc(left_arc);
  std::shared_ptr<DGArc> right_arc(new DGArcDirect(expr_vertex, right_oper));
  expr_vertex->add_exit_arc(right_arc);
#if DEBUG
  cout << "insert_expr_at: added arc between " << where->description()
       << " and " << expr_vertex->description() << endl;
#endif

  return expr_vertex;
}

// Replace recurrence relations with expressions
void DirectedGraph::remove_trivial_arithmetics() {
  using libint2::prefactor::Scalar;
  const std::shared_ptr<CTimeEntity<double> > const_one_point_zero =
      Scalar(1.0);
  const std::shared_ptr<CTimeEntity<double> > const_zero_point_zero =
      Scalar(0.0);
  typedef vertices::iterator iter;
  for (iter v = stack_.begin(); v != stack_.end(); ++v) {
    const ver_ptr& vptr = vertex_ptr(*v);
    std::shared_ptr<AlgebraicOperator<DGVertex> > oper_cast =
        std::dynamic_pointer_cast<AlgebraicOperator<DGVertex>, DGVertex>(
            (vptr));
    if (oper_cast) {
      // std::cout << "cast " << vptr->description() << " to " <<
      // oper_cast->description() << std::endl;
      typedef DGVertex::ArcSetType::const_iterator aciter;
      aciter a = oper_cast->first_exit_arc();
      std::shared_ptr<DGVertex> left = (*a)->dest();
      ++a;
      std::shared_ptr<DGVertex> right =
          oper_cast->num_exit_arcs() > 1
              ? (*a)->dest()
              : left;  // num_exit_arcs==1 is a corner case, e.g. 1*1

      using libint2::algebra::OperatorTypes;

      // 1.0 * x = x || 0.0 + x = x
      if (left->num_entry_arcs() == 1 &&
          ((oper_cast->type() == OperatorTypes::Times &&
            left->equiv(const_one_point_zero)) ||
           (oper_cast->type() == OperatorTypes::Plus &&
            left->equiv(const_zero_point_zero)))) {
#if DEBUG
        const bool success = remove_vertex_at((vptr), right);
        if (success) cout << "Removed vertex " << (vptr)->description() << endl;
#else
        remove_vertex_at((vptr), right);
#endif
      }
      // x * 1.0 = x || x + 0.0 = x
      else if (right->num_entry_arcs() == 1 &&
               ((oper_cast->type() == OperatorTypes::Times &&
                 right->equiv(const_one_point_zero)) ||
                (oper_cast->type() == OperatorTypes::Plus &&
                 right->equiv(const_zero_point_zero)))) {
#if DEBUG
        const bool success = remove_vertex_at((vptr), left);
        if (success) cout << "Removed vertex " << (vptr)->description() << endl;
#else
        remove_vertex_at((vptr), left);
#endif
      }

      // NOTE : more cases to come
    }
  }
}

namespace {
std::string stack_symbol(const std::shared_ptr<CodeContext>& ctext,
                         const DGVertex::Address& address,
                         const DGVertex::Size& size,
                         const std::string& low_rank, const std::string& veclen,
                         const std::string& prefix) {
  ostringstream oss;
  std::string stack_address = ctext->stack_address(address);
  oss << prefix << "[((hsi*" << size << "+" << stack_address << ")*" << low_rank
      << "+lsi)*" << veclen << "]";
  return oss.str();
}

/// Returns a "vector" form of stack symbol, e.g. converts libint->stack[x] to
/// libint->stack[x+vi]
inline std::string to_vector_symbol(const std::shared_ptr<DGVertex>& v) {
  std::string::size_type current_pos = 0;
  std::string symb = v->symbol();
  // replace repeatedly until the string is exhausted
  while (current_pos != std::string::npos) {
    // find "[" first
    const std::string left_braket("[");
    std::string::size_type where = symb.find(left_braket, current_pos);
    current_pos = where;
    // if the prefix indicating a stack symbol found:
    // 1) make sure vi doesn't appear between the brakets
    // 2) replace "]" with "+vi]"
    if (where != std::string::npos) {
      const std::string right_braket("]");
      std::string::size_type where = symb.find(right_braket, current_pos);
      if (where == std::string::npos)
        throw logic_error(
            "to_vector_symbol() -- address is set but no right braket found");

      const std::string forbidden("vi");
      std::string::size_type pos = symb.find(forbidden, current_pos);
      if (pos == std::string::npos || pos > where) {
        const std::string what_to_add("+vi");
        symb.insert(where, what_to_add);
        current_pos = where + 4;
      } else {
        current_pos = where + 1;
      }
    }
  }  // end of while
  return symb;
}
};  // namespace

//
// Handles "trivial" nodes. A node is trivial is it satisfies the following
// conditions: 0) not a target 1) has only one child 2) the exit arc is of a
// trivial type (DGArcDirect or IntegralSet_to_Integral applied to node of size
// 1)
//
// By "handling" I mean either removing the node from the graph or making a node
// refer to another node so that no code is generated for it.
//
void DirectedGraph::handle_trivial_nodes(
    const std::shared_ptr<CodeContext>& context) {
  typedef vertices::iterator iter;
  for (iter v = stack_.begin(); v != stack_.end(); ++v) {
    const ver_ptr& vptr = vertex_ptr(*v);
    // or if has more than 1 child
    if ((vptr)->num_exit_arcs() != 1) continue;
    std::shared_ptr<DGArc> arc = *((vptr)->first_exit_arc());

    // Is the exit arc DGArcRel<IntegralSet_to_Integrals> and (vptr)->size() ==
    // 1?
    if ((vptr)->size() == 1) {
      std::shared_ptr<DGArcRR> arc_cast =
          std::dynamic_pointer_cast<DGArcRR, DGArc>(arc);
      if (arc_cast) {
        std::shared_ptr<RecurrenceRelation> rr = arc_cast->rr();
        std::shared_ptr<IntegralSet_to_Integrals_base> rr_cast =
            std::dynamic_pointer_cast<IntegralSet_to_Integrals_base,
                                      RecurrenceRelation>(rr);
        if (rr_cast) {
          std::shared_ptr<DGVertex> child = arc->dest();

          // if (child->symbol_set() == false)
          {
            const std::string stack_name("stack");
            const std::shared_ptr<ImplicitDimensions>& dims =
                ImplicitDimensions::default_dims();
            std::string low_rank = dims->low_label();
            std::string veclen = dims->vecdim_label();

            if ((vptr)->address_set()) {
              child->set_symbol(stack_symbol(context, (vptr)->address() + 0,
                                             (vptr)->size(), low_rank, veclen,
                                             stack_name));
            }
            if ((vptr)->symbol_set()) {
              child->set_symbol(stack_symbol(context, 0, (vptr)->size(),
                                             low_rank, veclen,
                                             (vptr)->symbol()));
            }
          }

          vptr->refer_this_to(child);
          vptr->reset_symbol();
        }
      }
    }

    // Is the exit arc DGArcDirect?
    {
      std::shared_ptr<DGArcDirect> arc_cast =
          std::dynamic_pointer_cast<DGArcDirect, DGArc>(arc);
      if (arc_cast) {
        // remove the vertex, if possible
        // if this is a target -- cannot remove
        if (!(vptr)->is_a_target()) remove_vertex_at((vptr), arc->dest());
      }
    }

    // NOTE : more cases to come
  }
}

// If v1 and v2 are connected by DGArcDirect and all entry arcs to v1 are of the
// DGArcDirect type as well, this function will reattach all arcs extering v1 to
// v2 and remove v1 from the graph alltogether. return true if successful, false
// otherwise
bool DirectedGraph::remove_vertex_at(const std::shared_ptr<DGVertex>& v1,
                                     const std::shared_ptr<DGVertex>& v2) {
#if DEBUG
  cout << "remove_vertex_at: replacing " << v1->description() << " with "
       << v2->description() << endl;
#endif

  // Collect all entry arcs in a container
  DGVertex::ArcSetType v1_entry;
  typedef DGVertex::ArcSetType::iterator aiter;
  typedef DGVertex::ArcSetType::const_iterator aciter;
  const aciter abegin = v1->first_entry_arc();
  const aciter aend = v1->plast_entry_arc();
  // Verify that all entry arcs are DGArcDirect
  for (aciter a = abegin; a != aend; ++a) {
    // See if this is a direct arc -- otherwise cannot do this
    std::shared_ptr<DGArc> arc = (*a);
    std::shared_ptr<DGArcDirect> arc_cast =
        std::dynamic_pointer_cast<DGArcDirect, DGArc>(arc);
    if (arc_cast == 0) return false;
    v1_entry.push_back(*a);
#if DEBUG
    std::cout << "remove_vertex_at: examined v1 entry arc: from "
              << (*a)->orig()->description() << " to "
              << (*a)->dest()->description() << std::endl;
#endif
  }

  // Verify that v1 and v2 are connected by an arc and it is the only arc
  // exiting v1
  /*if (v1->num_exit_arcs() != 1 || v1->exit_arc(0)->dest() != v2)
    return false;*/

  // See if this is a direct arc -- otherwise cannot do this
  std::shared_ptr<DGArc> arc = *(v1->first_exit_arc());
  std::shared_ptr<DGArcDirect> arc_cast =
      std::dynamic_pointer_cast<DGArcDirect, DGArc>(arc);
  if (arc_cast == 0) return false;

    //
    // OK, now do work!
    //

#if DEBUG || DEBUG_RESTRUCTURE
  std::cout << "remove_vertex_at: v1" << endl;
  v1->print(std::cout);
  std::cout << "remove_vertex_at: v2" << endl;
  v2->print(std::cout);
#endif

  // Reconnect each of v1's entry arcs to v2
  unsigned int c = 0;
  for (aiter i = v1_entry.begin(); i != v1_entry.end(); ++i, ++c) {
    std::shared_ptr<DGVertex> parent = (*i)->orig();
#if DEBUG || DEBUG_RESTRUCTURE
    cout << "remove_vertex_at: replacing arc " << c << " connecting "
         << parent->description() << " to " << (*i)->dest()->description()
         << endl;
    cout << "remove_vertex_at: replacing arc " << c << " connecting " << parent
         << " to " << (*i)->dest() << endl;
#endif
    std::shared_ptr<DGArcDirect> new_arc(new DGArcDirect(parent, v2));
#if DEBUG || DEBUG_RESTRUCTURE
    cout << "remove_vertex_at:      with arc "
         << " connecting " << parent->description() << " to "
         << v2->description() << endl;
    cout << "remove_vertex_at:      with arc "
         << " connecting " << parent << " to " << v2 << endl;
#endif
    parent->replace_exit_arc(*i, new_arc);
#if DEBUG || DEBUG_RESTRUCTURE
    cout << "Replaced arcs: parent " << parent->description()
         << " now connected to " << new_arc->dest()->description() << endl;
    cout << "                ptr = " << parent << endl;
    const unsigned int nchildren = parent->num_exit_arcs();
    cout << "               parent has " << nchildren << " children" << endl;
    unsigned int c = 0;
    for (aciter a = parent->first_exit_arc(); a != parent->plast_exit_arc();
         ++a, ++c) {
      cout << "               child " << c << " " << (*a)->dest()->description()
           << endl;
      cout << "               child " << c << " ptr = " << (*a)->dest() << endl;
    }
#endif
  }

#if DEBUG || DEBUG_RESTRUCTURE
  std::cout << "remove_vertex_at: v1" << endl;
  v1->print(std::cout);
  std::cout << "remove_vertex_at: v2" << endl;
  v2->print(std::cout);
#endif

  // and fully disconnect this vertex
  v1->detach();
#if DEBUG || DEBUG_RESTRUCTURE
  std::cout << "remove_vertex_at: detached " << v1->description() << endl;
#endif

  return true;
}

void DirectedGraph::remove_disconnected_vertices() {
  typedef vertices::iterator iter;
  for (iter v = stack_.begin(); v != stack_.end();) {
    const ver_ptr& vptr = vertex_ptr(*v);
    iter vnext = v;
    ++vnext;  // note the next value of iterator before trying to erase
    if ((vptr)->num_entry_arcs() == 0 && (vptr)->num_exit_arcs() == 0 &&
        !(vptr)->is_a_target()) {
#if DEBUG
      cout << "Trying to erase disconnected vertex " << (vptr)->description()
           << " num_vertices = " << num_vertices() << endl;
#endif
      try {
        del_vertex(v);
      } catch (CannotPerformOperation& v) {
#if DEBUG
        cout << "But couldn't!!!" << endl;
#endif
        throw v;
      }
    }
    // update the iterator
    v = vnext;
  }
}

// generate_code uses this helper function. It's declared in libint2 namespace
// because it's also used elsewhere
namespace libint2 {
std::string declare_function(const std::shared_ptr<CodeContext>& context,
                             const std::shared_ptr<ImplicitDimensions>& dims,
                             const std::shared_ptr<CodeSymbols>& args,
                             const std::string& tlabel,
                             const std::string& function_descr,
                             std::ostream& decl) {
  std::string function_name = label_to_funcname(function_descr);
  function_name = context->label_to_name(function_name);

  decl << context->code_prefix();
  std::string func_decl;
  std::ostringstream oss;
  oss << context->type_name<void>() << " " << function_name << "("
      << context->const_modifier() << context->inteval_type_name(tlabel)
      << "* inteval";
  const unsigned int nargs = args->n();
  if (nargs > 0) {
    // first argument is always the target, which is never a const
    oss << ", " << context->type_name<double*>() << " " << args->symbol(0);
    for (unsigned int a = 1; a < nargs; a++) {
      oss << ", " << context->type_name<const double*>() << " "
          << args->symbol(a);
    }
  }
  if (!dims->high_is_static()) {
    oss << ", " << context->type_name<int>() << " " << dims->high()->id();
  }
  if (!dims->low_is_static()) {
    oss << ", " << context->type_name<int>() << " " << dims->low()->id();
  }
  oss << ")";
  func_decl = oss.str();

  decl << func_decl << context->end_of_stat() << endl;
  decl << context->code_postfix();

  return func_decl;
}
}  // namespace libint2

namespace {
template <class Container>
std::shared_ptr<MemBlockSet> to_memoryblks(Container& vertices) {
  std::shared_ptr<MemBlockSet> result(new MemBlockSet);
  typedef typename Container::const_iterator citer;
  typedef typename Container::iterator iter;
  citer end(vertices.end());
  for (iter v = vertices.begin(); v != end; ++v) {
    const DirectedGraph::ver_ptr& vptr = vertex_ptr(*v);
    result->push_back(MemBlock((vptr)->address(), (vptr)->size(), false,
                               std::shared_ptr<MemBlock>(),
                               std::shared_ptr<MemBlock>()));
  }
  return result;
}
};  // namespace

//
//
//
void DirectedGraph::generate_code(
    const std::shared_ptr<CodeContext>& context,
    const std::shared_ptr<MemoryManager>& memman,
    const std::shared_ptr<ImplicitDimensions>& dims,
    const std::shared_ptr<CodeSymbols>& args, const std::string& label,
    std::ostream& decl, std::ostream& def) {
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  const std::string tlabel = taskmgr.current().label();

  decl << context->copyright();
  decl << context->std_header();
  std::string comment("This code computes ");
  comment += label;
  comment += "\n";
  if (context->comments_on()) decl << context->comment(comment) << endl;

  const std::string func_decl =
      declare_function(context, dims, args, tlabel, label, decl);

  // are there prerequisite vertices that are not precomputed? Then may need to
  // call precompute function
  const bool missing_prereqs = this->missing_prerequisites();
  std::string func_prereq_name;
  std::string func_prereq_decl;
  if (missing_prereqs) {
    std::ostringstream oss;
    func_prereq_name =
        context->label_to_name(label_to_funcname(label)) + "_prereq";
    func_prereq_decl =
        declare_function(context, dims, args, tlabel, func_prereq_name, oss);
  }

  //
  // Generate function's definition
  //

  def << context->copyright();
  // include standard headers
  def << context->std_header();
  // include declarations for all function calls:
  // 1) update func_names_
  // 2) (optional) if will compute prerequisites add the name of the function
  // that will compute them 3) include their headers into the current definition
  // file
  update_func_names();
  if (missing_prereqs) func_names_[func_prereq_name] = true;
  for (FuncNameContainer::const_iterator fn = func_names_.begin();
       fn != func_names_.end(); fn++) {
    string function_name = (*fn).first;
    def << "#include <"
        << context->label_to_name(context->cparams()->api_prefix() +
                                  function_name)
        << ".h>" << endl;
  }
  def << endl;

  def << context->code_prefix();
  def << func_decl << context->open_block() << endl;
  def << context->std_function_header();

  // allocate data and assign symbols
  context->reset();
  // if we vectorize by-line then all data is allocated on Libint's stack
  if (context->cparams()->vectorize_by_line()) allocate_mem(memman, dims, 0);
  // otherwise only arrays go on Libint's stack (scalars are handled by the
  // compiler)
  else
    allocate_mem(memman, dims, 1);
  assign_symbols(context, dims);

  // then ...
  if (missing_prereqs) {  // need to compute prerequisites?

    // if profiling is on, start the next timer, after stopping the current
    // timer (if possible)
    if (registry()->current_timer() >= 0) {
      const int current_timer = registry()->current_timer();
      def << context->macro_if("LIBINT2_CPLUSPLUS_STD >= 2011");
      if (current_timer > 0)
        def << "inteval->timers->stop(" << current_timer << ");" << std::endl;
      def << "inteval->timers->start(" << current_timer + 1 << ");"
          << std::endl;
      def << context->macro_endif();
    }

    //
    // need to zero out space for all missing prerequisites -- prereq evaluator
    // will accumulate into that space
    //

    // get prerequisites
    PrerequisitesExtractor pe;
    this->foreach (pe);
    // merge their blocks (likely into one -- allocate_mem should have taken
    // care of that)
    std::shared_ptr<MemBlockSet> targetblks = to_memoryblks(pe.vertices);
    merge(*targetblks);

    // zero each one out
    for (MemBlockSet::iterator b = targetblks->begin(); b != targetblks->end();
         ++b) {
      const size s = b->size();
      std::shared_ptr<CTimeEntity<int> > bdim(new CTimeEntity<int>(s));

      std::shared_ptr<Entity> bvecdim;
      if (!dims->vecdim_is_static()) {
        std::shared_ptr<RTimeEntity<EntityTypes::Int> > vecdim =
            std::dynamic_pointer_cast<RTimeEntity<EntityTypes::Int>, Entity>(
                dims->vecdim());
        bvecdim = vecdim * bdim;
      } else {
        std::shared_ptr<CTimeEntity<int> > vecdim =
            std::dynamic_pointer_cast<CTimeEntity<int>, Entity>(dims->vecdim());
        bvecdim = vecdim * bdim;
      }
      def << "_libint2_static_api_bzero_short_(" << registry()->stack_name()
          << "+" << b->address() << "*" << dims->vecdim()->id() << ","
          << bvecdim->id() << ")" << endl;
    }

    // and call the prereq evaluator
    std::shared_ptr<Entity> zero(new CTimeEntity<int>(0));
    std::shared_ptr<Entity> contr_depth(
        new RTimeEntity<EntityTypes::Int>("contrdepth"));
    std::string contr_index("c");
    std::shared_ptr<ForLoop> contr_loop(
        new ForLoop(context, contr_index, contr_depth, zero));

    def << context->decldef("const int", "contrdepth", "inteval->contrdepth");
    def << contr_loop->open();
    def << func_prereq_name << "(inteval+c, " << registry()->stack_name() << ")"
        << context->end_of_stat() << endl;
    def << contr_loop->close() << endl;

    // if profiling is on, start the next timer, after stopping the current
    // timer (if possible)
    if (registry()->current_timer() >= 0) {
      const int current_timer = registry()->current_timer();
      def << context->macro_if("LIBINT2_CPLUSPLUS_STD >= 2011");
      def << "inteval->timers->stop(" << current_timer + 1 << ");" << std::endl;
      if (current_timer > 0)
        def << "inteval->timers->start(" << current_timer << ");" << std::endl;
      def << context->macro_endif();
    }
  }

  // if profiling=on and the current timer is outermost, start it
  // if this is not outermost timer, it has been started outside of this
  // function
  if (registry()->current_timer() == 0) {
    def << context->macro_if("LIBINT2_CPLUSPLUS_STD >= 2011");
    def << "inteval->timers->start(0);" << std::endl;
    def << context->macro_endif();
  }

  // now print out the code for this graph
  print_def(context, def, dims, args);

  // if profiling=on and the current timer is outermost, stop it
  // if this is not outermost timer, it is managed outside of this function
  if (registry()->current_timer() == 0) {
    def << context->macro_if("LIBINT2_CPLUSPLUS_STD >= 2011");
    def << "inteval->timers->stop(0);" << std::endl;
    def << context->macro_endif();
  }

  def << context->close_block() << endl;
  def << context->code_postfix();
}

void DirectedGraph::allocate_mem(
    const std::shared_ptr<MemoryManager>& memman,
    const std::shared_ptr<ImplicitDimensions>& dims,
    unsigned int min_size_to_alloc) {
  // NOTE does this belong here?
  // First, reset tag counters
  prepare_to_traverse();

  struct TargetAllocator {
    typedef DirectedGraph::targets::const_iterator target_citer;
    typedef DirectedGraph::size sz;

    const DirectedGraph::targets& targets_;
    const std::shared_ptr<MemoryManager>& memman_;
    bool all_targets_;
    sz size_;

    TargetAllocator(const DirectedGraph::targets& t,
                    const std::shared_ptr<MemoryManager>& mm, bool all_targets)
        : targets_(t), memman_(mm), all_targets_(all_targets) {
      // compute the aggregate size of all targets
      target_citer end = targets_.end();
      size_ = 0;
      for (target_citer t = targets_.begin(); t != end; ++t) {
        const ver_ptr& tptr = vertex_ptr(*t);
        if (all_targets_ || (!tptr->symbol_set() && !tptr->address_set())) {
          size_ += (tptr)->size();
        }
      }
    }

    sz size() const { return size_; }

    void allocate() {
      for (target_citer v = targets_.begin(); v != targets_.end(); ++v) {
        const ver_ptr& vptr = vertex_ptr(*v);
        if (all_targets_ || (!vptr->symbol_set() && !vptr->address_set())) {
          const MemoryManager::Address address = memman_->alloc(vptr->size());
          // if the vertex is just an alias, pass the address on
          if (vptr->refers_to_another())
            (*(vptr->first_exit_arc()))->dest()->set_address(address);
          else
            vptr->set_address(address);
        }
      }
    }
  };

  //
  // First, allocate all prerequisites that are not precomputed
  // since they will be computed before the evaluation of this graph
  // they will be targets of the previous computation and will be
  // at the beginning of the previous stack. Preallocate them here.
  //
  if (this->missing_prerequisites()) {
    PrerequisitesExtractor pe;
    this->foreach (pe);
    targets prereqs(pe.vertices.size());
    std::copy(pe.vertices.begin(), pe.vertices.end(), prereqs.begin());
    // prereqs will be put on the prereq DirectedGraph in the same order
    // so use the same allocation mechanism here as for targets
    const bool all_targets = true;
    TargetAllocator ta(prereqs, memman, all_targets);
    ta.allocate();
  }

  //
  // If need to accumulate targets, special events must happen here.
  //
  // NOTES on how to handle accumulation
  // 1) if all targets are unrolled then need to identify which integrals are
  // part of target sets and use += instead of = 2) if no targets are unrolled
  // then allocate extra space for the target quartets and, after code has been
  // generated,
  //    accumulate target sets into those
  //    EXCEPTION non new space for integral sets is needed if they are
  //    decontracted
  // 3) if some targets are not unrolled then still need the extra space. The
  // targets which were not unrolled should be handled
  //    as usual, i.e. not accumulated -- accumulation happens at the end
  if (registry()->accumulate_targets()) {
    // need extra buffers for targets if some are not unrolled ( and not
    // decontracted)
    // const bool need_copies_of_targets = nonunrolled_targets(targets_);
    const bool need_copies_of_targets =
        std::find_if(
            targets_.begin(), targets_.end(), [](std::shared_ptr<DGVertex> i) {
              return !DecontractedIntegralSet()(i) && !UnrolledIntegralSet()(i);
            }) != targets_.end();

    iregistry()->accumulate_targets_directly(!need_copies_of_targets);

    if (need_copies_of_targets) {
      const bool all_targets = true;
      TargetAllocator ta(targets_, memman, all_targets);
      const size size_of_targets = ta.size();
      iregistry()->size_of_target_accum(size_of_targets);

      // allocate every target accumulator manually
      const address targets_buffer = memman->alloc(size_of_targets);
      address curr_ptr = targets_buffer;
      for (target_citer t = targets_.begin(); t != targets_.end(); ++t) {
        const ver_ptr& tptr = vertex_ptr(*t);
        target_accums_.push_back(curr_ptr);
        curr_ptr += (tptr)->size();
      }

    }  // need copies of targets
  }    // need to accumulate targets

  // Second, MUST allocate space for all targets whose symbols are not set
  // explicitly If a symbol is set means the object is not on stack (e.g. if
  // location of target is passed as an argument to set-level function) This
  // code ensures that target quartets are persistent, i.e. never overwritten,
  // and can be accumulated into
  {
    const bool all_targets = false;  // only targets without symbols
    TargetAllocator ta(targets_, memman, all_targets);
    ta.allocate();
  }

  //
  // How memory management happens:
  // Go through the traversal order and at each step tag every child
  // Once a child receives same number of tags as the number of parents,
  // it can be deallocated
  //
  std::shared_ptr<DGVertex> vertex = first_to_compute_;
  do {
    std::shared_ptr<DGArcRR> arcrr;
    // memory only needs to be managed for some quantities:
    // this conditional decides whether this vertex is on the stack
    bool need_to_allocate = true;

    // If symbol is set then the object is not on stack
    need_to_allocate &= not vertex->symbol_set();

    // if address is already set, no need to manage
    need_to_allocate &= not vertex->address_set();

    // precomputed objects don't go on stack
    need_to_allocate &= not vertex->precomputed();

    // don't allocate on stack if smaller than or equal to min_size_to_alloc
    // two exceptions, however:
    // 1) it's a target
    // 2) it's an unrolled integral set whose members are not precomputed
    // quantities
    //    typically integral sets of size 1 are precomputed and don't need to be
    //    on stack, however if they are not, the integral will need to be stored
    //    somewhere and the rule for assigning code symbols to members of
    //    unrolled integral sets requires the integral set to have an address
    //    assigned
    need_to_allocate &=
        (vertex->size() > min_size_to_alloc || vertex->is_a_target() ||
         (vertex->size() == vertex->num_exit_arcs() &&
          ((arcrr = std::dynamic_pointer_cast<DGArcRR, DGArc>(
                *(vertex->first_exit_arc()))) != 0
               ? std::dynamic_pointer_cast<IntegralSet_to_Integrals_base,
                                           RecurrenceRelation>(arcrr->rr()) != 0
               : false) &&
          !(*(vertex->first_exit_arc()))->dest()->precomputed()));

    // if this is a member of unrolled integral set whose address or symbol is
    // set, no need to allocate
    if (need_to_allocate && vertex->num_entry_arcs() != 0) {
      arcrr = std::dynamic_pointer_cast<DGArcRR, DGArc>(
          *(vertex->first_entry_arc()));
      if (arcrr) {
        if (std::dynamic_pointer_cast<IntegralSet_to_Integrals_base,
                                      RecurrenceRelation>(arcrr->rr()) != 0) {
          if (arcrr->orig()->symbol_set() || arcrr->orig()->address_set())
            need_to_allocate = false;
        }
      }
    }

    if (need_to_allocate) {
      MemoryManager::Address addr = memman->alloc(vertex->size());
      vertex->set_address(addr);
#if DEBUG
      cout << "allocated " << vertex->description() << " at " << addr
           << " (size=" << vertex->size() << ")" << endl;
#endif

      typedef DGVertex::ArcSetType::const_iterator aciter;
      const aciter abegin = vertex->first_exit_arc();
      const aciter aend = vertex->plast_exit_arc();
      // Verify that all entry arcs are DGArcDirect
      for (aciter a = abegin; a != aend; ++a) {
        std::shared_ptr<DGVertex> child = (*a)->dest();
        const unsigned int ntags = child->tag();
        // Do NOT deallocate if it's a target!
        if (ntags == child->num_entry_arcs() && child->address_set() &&
            !child->is_a_target()) {
          memman->free(child->address());
        }
      }
    }
    vertex = vertex->postcalc();
  } while (vertex != 0);
}

void DirectedGraph::assign_symbols(
    const std::shared_ptr<CodeContext>& context,
    const std::shared_ptr<ImplicitDimensions>& dims) {
  std::ostringstream os;
  const std::string null_str("");
  const std::string stack_name = registry()->stack_name();
  // There used to be a compiler/library/bug on OS X? The above failed...
  // Valgrind under Linux shows no memory problems...
  // const std::string stack_name("stack");

  // Generate the label for the rank of the low dimension
  std::string low_rank = dims->low_label();
  std::string veclen = dims->vecdim_label();

  // First, set symbols for all vertices which have address assigned
  typedef vertices::iterator iter;
  for (iter v = stack_.begin(); v != stack_.end(); ++v) {
    const ver_ptr& vptr = vertex_ptr(*v);
    if (!(vptr)->symbol_set() && (vptr)->address_set()) {
      (vptr)->set_symbol(stack_symbol(context, (vptr)->address(),
                                      (vptr)->size(), low_rank, veclen,
                                      stack_name));
    }
  }

  // Second, find all nodes which were unrolled using IntegralSet_to_Integrals:
  // 1) such nodes do not need symbols generated since they never appear in the
  // code expicitly 2) children of such nodes have symbols that depend on the
  // parent's address upstream such nodes of size one were aliased ("referred")
  // to their children -- just skip these
  for (iter v = stack_.begin(); v != stack_.end(); ++v) {
    const ver_ptr& vptr = vertex_ptr(*v);
    if ((vptr)->num_exit_arcs() == 0) continue;
    if (vptr->refers_to_another()) continue;
    std::shared_ptr<DGArc> arc = *((vptr)->first_exit_arc());
    std::shared_ptr<DGArcRR> arc_rr =
        std::dynamic_pointer_cast<DGArcRR, DGArc>(arc);
    if (arc_rr == 0) continue;
    std::shared_ptr<RecurrenceRelation> rr = arc_rr->rr();
    std::shared_ptr<IntegralSet_to_Integrals_base> iset_to_i =
        std::dynamic_pointer_cast<IntegralSet_to_Integrals_base,
                                  RecurrenceRelation>(rr);
    if (iset_to_i == 0) {
      continue;
    } else {
      typedef DGVertex::ArcSetType::const_iterator aciter;
      const aciter abegin = (vptr)->first_exit_arc();
      const aciter aend = (vptr)->plast_exit_arc();
      unsigned int c = 0;
      // Verify that all entry arcs are DGArcDirect
      for (aciter a = abegin; a != aend; ++a, ++c) {
        std::shared_ptr<DGVertex> child = (*a)->dest();
        // If the child is precomputed and it's parent symbol is not set -- its
        // symbol will be set as usual
        if (!child->precomputed() || (vptr)->symbol_set()) {
          // check if symbol is already set ... this indicates interference with
          // the logic upstream, it should be set here?
          if (child->symbol_set()) {
            std::cout << "WARNING: symbol for " << child->description()
                      << " (unrolled integral) already set" << std::endl;
            continue;
          }

          if ((vptr)->address_set()) {
            child->set_symbol(stack_symbol(context, (vptr)->address() + c,
                                           (vptr)->size(), low_rank, veclen,
                                           stack_name));
          } else {
            child->set_symbol(stack_symbol(context, c, (vptr)->size(), low_rank,
                                           veclen, (vptr)->symbol()));
          }
        }
      }
      (vptr)->refer_this_to((*((vptr)->first_exit_arc()))->dest());
      (vptr)->reset_symbol();
    }
  }

  // then process all other symbols, EXCEPT operators
  for (iter v = stack_.begin(); v != stack_.end(); ++v) {
    const ver_ptr& vptr = vertex_ptr(*v);
#if DEBUG
    cout << "Trying to assign symbol to " << (vptr)->description() << endl;
#endif
    if ((vptr)->symbol_set()) {
#if DEBUG
      cout << "symbol already set to " << (vptr)->symbol() << endl;
#endif
      continue;
    }

    // test if the vertex is a static quantity, like a constant
    {
      typedef CTimeEntity<double> cdouble;
      std::shared_ptr<cdouble> ptr_cast =
          std::dynamic_pointer_cast<cdouble, DGVertex>((vptr));
      if (ptr_cast) {
        (vptr)->set_symbol(ptr_cast->label());
        continue;
      }
    }

    // test if the vertex is precomputed runtime quantity, like a geometric
    // parameter
    if ((vptr)->precomputed()) {
      std::string symbol("inteval->");
      symbol += context->label_to_name((vptr)->label());
      symbol += "[vi]";
      (vptr)->set_symbol(symbol);
      continue;
    }

    // test if the vertex is other runtime quantity
    {
      typedef RTimeEntity<double> cdouble;
      std::shared_ptr<cdouble> ptr_cast =
          std::dynamic_pointer_cast<cdouble, DGVertex>((vptr));
      if (ptr_cast) {
        (vptr)->set_symbol(ptr_cast->label());
        continue;
      }
    }

    // test if the vertex is some other data that was not allocated on stack
#if 1
    if (vptr->size() == 1) {
      {  // basic integral
        typedef IntegralSet<IncableBFSet> intset;
        std::shared_ptr<intset> ptr_cast =
            std::dynamic_pointer_cast<intset, DGVertex>((vptr));
        if (ptr_cast) {
          (vptr)->set_symbol(context->unique_name<EntityTypes::FP>());
          continue;
        }
      }
    }
#endif

  }  // done with everything BUT operators

  // finally, process all operators (start with most recently added vertices
  // since those are much more likely to be on the bottom of the graph).
  typedef vertices::reverse_iterator riter;
  for (riter v = stack_.rbegin(); v != stack_.rend(); ++v) {
    ver_ptr& vptr = vertex_ptr(*v);
    if (vptr->symbol_set()) continue;
#if DEBUG
    cout << "Trying to assign symbol to operator " << (vptr)->description()
         << endl;
#endif
    assign_oper_symbol(context, (vptr));
  }
}

void DirectedGraph::assign_oper_symbol(
    const std::shared_ptr<CodeContext>& context,
    std::shared_ptr<DGVertex>& vertex) {
  // do nothing if the vertex has a symbol or is not an operator
  if (vertex->symbol_set()) return;

  {
    typedef AlgebraicOperator<DGVertex> oper;
    std::shared_ptr<oper> ptr_cast =
        std::dynamic_pointer_cast<oper, DGVertex>(vertex);
    if (ptr_cast) {
      // is it in a subtree?
      const bool on_a_subtree = (vertex->subtree() != 0);

      // If no -- it will be an automatic variable
      if (!on_a_subtree)
        vertex->set_symbol(context->unique_name<EntityTypes::FP>());
      // else assign symbols to left and right arguments
      else {
        typedef DGVertex::ArcSetType::const_iterator aciter;
        aciter arc = ptr_cast->first_exit_arc();
        std::shared_ptr<DGVertex> left = (*arc)->dest();
        ++arc;
        std::shared_ptr<DGVertex> right = (*arc)->dest();
        assign_oper_symbol(context, left);
        assign_oper_symbol(context, right);

        std::ostringstream oss;
        oss << "( " << left->symbol() << " ) " << ptr_cast->label() << " ( "
            << right->symbol() << " )";
        vertex->set_symbol(oss.str());
      }
    }
  }
}

namespace {
/// returns how many times token appears in str
unsigned int nfind(const std::string& str, const std::string& token) {
  unsigned int nfinds = 0;
  typedef std::string::size_type size_type;
  size_type current_pos = 0;
  while (1) {
    current_pos = str.find(token, current_pos);
    if (current_pos != std::string::npos) {
      ++nfinds;
      ++current_pos;
    } else
      return nfinds;
  }
}

#define DO_NOT_COUNT_DIV 1
/// Returns the number of FLOPs in an expression
unsigned int nflops(const std::string& expr) {
  static const std::string mul(" * ");
  static const std::string div(" / ");
  static const std::string plus(" + ");
  static const std::string minus(" - ");
  const unsigned int nflops = nfind(expr, mul) + nfind(expr, plus) +
                              nfind(expr, minus)
#if !DO_NOT_COUNT_DIV
                              + nfind(expr, div)
#endif
      ;
  return nflops;
}
}  // namespace

void DirectedGraph::print_def(const std::shared_ptr<CodeContext>& context,
                              std::ostream& os,
                              const std::shared_ptr<ImplicitDimensions>& dims,
                              const std::shared_ptr<CodeSymbols>& args) {
  std::ostringstream oss;
  const std::string null_str("");
  std::shared_ptr<Entity> ctimeconst_zero(new CTimeEntity<int>(0));

  //
  // set stack ... if this function was given any arguments, the first is always
  // stack
  //
  os << context->decldef(
      context->type_name<double* const>(), "stack",
      (args->n() >= 1) ? args->symbol(0) : registry()->stack_name());

  const bool accumulate_targets_directly =
      registry()->accumulate_targets() &&
      iregistry()->accumulate_targets_directly();
  const bool accumulate_targets_indirectly =
      registry()->accumulate_targets() && !accumulate_targets_directly;

  //
  // To optimize accumulation and setting blocks to zero, need to merge
  // maximally the targets (the accumulation area is one block)
  //
  std::shared_ptr<MemBlockSet> targetblks;
  if (registry()->accumulate_targets()) {
    if (accumulate_targets_directly) {
      targetblks = to_memoryblks(targets_);
      merge(*targetblks);
    } else {
      targetblks = std::shared_ptr<MemBlockSet>(new MemBlockSet);
      targetblks->push_back(MemBlock(0, iregistry()->size_of_target_accum(),
                                     false, std::shared_ptr<MemBlock>(),
                                     std::shared_ptr<MemBlock>()));
    }
  }

  //
  // If accumulating integrals, check inteval's zero_out_targets. If set to 1 --
  // zero out accumulated targets
  //
#if LIBINT_ACCUM_INTS
  if (registry()->accumulate_targets()) {
    const bool vecdim_is_static = dims->vecdim_is_static();

    os << "if (inteval->zero_out_targets) {" << std::endl;

    typedef MemBlockSet::const_iterator citer;
    typedef MemBlockSet::iterator iter;
    const citer end = targetblks->end();
    for (iter b = targetblks->begin(); b != end; ++b) {
      size s = b->size();
      std::shared_ptr<CTimeEntity<int> > bdim(new CTimeEntity<int>(s));

      std::shared_ptr<Entity> bvecdim;
      if (!vecdim_is_static) {
        std::shared_ptr<RTimeEntity<EntityTypes::Int> > vecdim =
            std::dynamic_pointer_cast<RTimeEntity<EntityTypes::Int>, Entity>(
                dims->vecdim());
        bvecdim = vecdim * bdim;
      } else {
        std::shared_ptr<CTimeEntity<int> > vecdim =
            std::dynamic_pointer_cast<CTimeEntity<int>, Entity>(dims->vecdim());
        bvecdim = vecdim * bdim;
      }

#if 0
      std::string loopvar("i");
      ForLoop loop(context,loopvar,bvecdim,ctimeconst_zero);
      os << loop.open();
      {
        ostringstream oss;
        oss << registry()->stack_name() << "[" << loopvar << "]";
        const std::string zero("0");
        os << context->assign(oss.str(),zero);
      }
      os << loop.close();
#endif
      os << "_libint2_static_api_bzero_short_(" << registry()->stack_name()
         << "+" << b->address() << "*" << dims->vecdim()->id() << ","
         << bvecdim->id() << ")" << endl;
    }

    os << "inteval->zero_out_targets = 0;" << std::endl << "}" << std::endl;
  }
#endif

  std::string varname("hsi");
  std::shared_ptr<ForLoop> hsi_loop(
      new ForLoop(context, varname, dims->high(),
                  std::shared_ptr<Entity>(new CTimeEntity<int>(0))));
  os << hsi_loop->open();

  varname = "lsi";
  std::shared_ptr<ForLoop> lsi_loop(
      new ForLoop(context, varname, dims->low(),
                  std::shared_ptr<Entity>(new CTimeEntity<int>(0))));
  os << lsi_loop->open();

  // the vector loop is created outside of the body of the function if
  // 1) blockwise vectorization is requested
  // and
  // 2) this is a purely int-unit code, i.e. there are no explicit RRs on sets
  // in the body Otherwise, I create a dummy vector loop with the vector loop
  // index set to 0
  const unsigned int max_vector_length =
      context->cparams()->max_vector_length();
  const bool vectorize = (max_vector_length != 1);
  const bool vectorize_by_line = context->cparams()->vectorize_by_line();
  const bool create_outer_vector_loop =
      !vectorize_by_line && !cannot_enclose_in_outer_vloop();
  varname = "vi";
  // outer vector loop
  std::shared_ptr<ForLoop> outer_vloop;
  // vector loop for each code line
  std::shared_ptr<ForLoop> line_vloop;
  if (create_outer_vector_loop) {
    std::shared_ptr<ForLoop> tmp_vi_loop(
        new ForLoop(context, varname, dims->vecdim(),
                    std::shared_ptr<Entity>(new CTimeEntity<int>(0))));
    outer_vloop = tmp_vi_loop;
  } else {
    std::shared_ptr<Entity> unit_dim(new CTimeEntity<int>(1));
    std::shared_ptr<ForLoop> tmp_vi_loop(
        new ForLoop(context, varname, unit_dim,
                    std::shared_ptr<Entity>(new CTimeEntity<int>(0))));
    outer_vloop = tmp_vi_loop;
    std::shared_ptr<ForLoop> tmp2_vi_loop(
        new ForLoop(context, varname, dims->vecdim(),
                    std::shared_ptr<Entity>(new CTimeEntity<int>(0))));
    line_vloop = tmp2_vi_loop;
    // note that both loops use same variable name -- standard C++ scoping rules
    // allow it -- but the outer loop will become a declaration of a constant
    // variable
  }
  os << outer_vloop->open();

  //
  // generate code for vertices
  //
  unsigned int nflops_total = 0;
  std::shared_ptr<DGVertex> current_vertex = first_to_compute_;
  do {
    // skip if already scheduled
    if (current_vertex->scheduled()) goto next;

    // skip if this is a dummy vertex (i.e. refers to someone else)
    if (current_vertex->refers_to_another()) goto next;

    // for every vertex that has a defined symbol, hence must be defined in code
    if (current_vertex->symbol_set()) {
      // Type declaration if this is an automatic variable
      // ids for automatic variables cannot have characters '[' or ']'
      const std::string& symbol = current_vertex->symbol();
      if (symbol.find('[') == std::string::npos) {
        os << context->declare(context->type_name<double>(),
                               current_vertex->symbol());
#if CHECK_SAFETY
        current_vertex->declared(true);
#endif
      }

      // print algebraic expression
      {
        typedef AlgebraicOperator<DGVertex> oper_type;
        std::shared_ptr<oper_type> oper_ptr =
            std::dynamic_pointer_cast<oper_type, DGVertex>(current_vertex);
        if (oper_ptr) {
          // If this is an Integral in a target IntegralSet AND
          // can accumulate targets directly -- use '+=' instead of '='
          const bool accumulate_not_assign =
              accumulate_targets_directly &&
              IntegralInTargetIntegralSet()(current_vertex);

          typedef DGVertex::ArcSetType::const_iterator aciter;
          aciter a = oper_ptr->first_exit_arc();
          auto left_arg = (*a)->dest();
          ++a;
          auto right_arg = (*a)->dest();

          if (context->comments_on()) {
            oss.str(null_str);
            oss << current_vertex->label()
                << (accumulate_not_assign ? " += " : " = ") << left_arg->label()
                << oper_ptr->label() << right_arg->label();
            os << context->comment(oss.str()) << endl;
          }

          // expression
#if DEBUG
          cout << "Generating code for " << current_vertex->description()
               << endl;
          cout << "              ptr = " << current_vertex << endl;
          cout << "         left_arg = " << left_arg->description() << endl;
          cout << "        right_arg = " << right_arg->description() << endl;
#endif

          // can we generate a composite instruction, like FMA?
          // - FMA: yes if:
          //     o current_vertex is a multiply
          //     o it has one parent and the parent is a +/- operator
          //     o parent's other argument has been computed already
          bool generate_fma = false;
          std::shared_ptr<oper_type> parent_oper_ptr;
          std::shared_ptr<DGVertex> fma_other_arg;
#if LIBINT_GENERATE_FMA
          {
            if (oper_ptr->type() == algebra::OperatorTypes::Times &&
                oper_ptr->num_entry_arcs() == 1) {
              parent_oper_ptr = std::dynamic_pointer_cast<oper_type, DGVertex>(
                  (*(oper_ptr->first_entry_arc()))->orig());
              if (parent_oper_ptr != 0) {
                auto parent_oper_type = parent_oper_ptr->type();
                if (parent_oper_type == algebra::OperatorTypes::Plus ||
                    parent_oper_type == algebra::OperatorTypes::Minus) {
                  auto arg1 = (*(parent_oper_ptr->first_exit_arc()))->dest();
                  auto arg2 =
                      (*(++(parent_oper_ptr->first_exit_arc())))->dest();
                  auto other_arg = (arg1 == current_vertex) ? arg2 : arg1;
                  const bool other_ready =
                      other_arg->num_exit_arcs() == 0 || other_arg->scheduled();

                  // can do fmadd if other_ready
                  // can do fmsub if other_ready and it's arg2
                  if (parent_oper_type == algebra::OperatorTypes::Plus)
                    generate_fma = other_ready;
                  else
                    generate_fma = other_ready && (current_vertex == arg1);
                  if (generate_fma) {
                    // std::cout << context->comment("CAN GENERATE FMA!!!") <<
                    // endl;
                    fma_other_arg = other_arg;

#if DEBUG
                    std::cout << "Generating FMA:" << std::endl;
                    std::cout << "parent:\n";
                    parent_oper_ptr->print(std::cout);
                    std::cout << "current_vertex:\n";
                    current_vertex->print(std::cout);
                    std::cout << "other_vertex:\n";
                    other_arg->print(std::cout);
                    std::cout << "arg1:\n";
                    arg1->print(std::cout);
                    std::cout << "arg2:\n";
                    arg2->print(std::cout);
                    std::cout << "child1:\n";
                    (*(parent_oper_ptr->first_exit_arc()))
                        ->dest()
                        ->print(std::cout);
                    std::cout << "child2:\n";
                    (*(++(parent_oper_ptr->first_exit_arc())))
                        ->dest()
                        ->print(std::cout);
#endif

                    // may need to declare the variable for the parent
                    if (parent_oper_ptr->symbol_set()) {
                      // Type declaration if this is an automatic variable
                      // ids for automatic variables cannot have characters '['
                      // or ']'
                      const std::string& symbol = parent_oper_ptr->symbol();
                      if (symbol.find('[') == std::string::npos) {
                        os << context->declare(context->type_name<double>(),
                                               parent_oper_ptr->symbol());
#if CHECK_SAFETY
                        parent_oper_ptr->declared(true);
#endif
                      }
                    }
                  }
                }
              }
            }
          }
#endif  // LIBINT_GENERATE_FMA=1

          // convert symbols to their vector form if needed
          std::string curr_symbol = current_vertex->symbol();
          std::string left_symbol = left_arg->symbol();
          std::string right_symbol = right_arg->symbol();
          std::string parent_symbol =
              generate_fma ? parent_oper_ptr->symbol() : "";
          std::string fma_other_arg_symbol =
              generate_fma ? fma_other_arg->symbol() : "";
          if (vectorize) {
            curr_symbol = to_vector_symbol(current_vertex);
            left_symbol = to_vector_symbol(left_arg);
            right_symbol = to_vector_symbol(right_arg);
            if (generate_fma) {
              parent_symbol = to_vector_symbol(parent_oper_ptr);
              fma_other_arg_symbol = to_vector_symbol(fma_other_arg);
            }
          }
#if CHECK_SAFETY && 0
          bool left_not_declared =
              left_arg->need_to_compute() && !left_arg->declared();
          bool right_not_declared =
              right_arg->need_to_compute() && !right_arg->declared();
          if (left_not_declared || right_not_declared) {
            std::cout << "Current vertex:" << endl;
            current_vertex->print(std::cout);
            std::cout << "left arg      :" << endl;
            left_arg->print(std::cout);
            std::cout << "right arg     :" << endl;
            right_arg->print(std::cout);
            if (left_not_declared)
              throw ProgrammingError(
                  "DirectedGraph::print_def() -- left_arg not declared");
            if (right_not_declared)
              throw ProgrammingError(
                  "DirectedGraph::print_def() -- right_arg not declared");
          }
#endif

          if (vectorize_by_line) os << line_vloop->open();
          // the statement that does the work
          {
            if (accumulate_not_assign) {
              if (generate_fma) {
                os << context->accumulate_ternary_expr(
                    parent_symbol, left_symbol, oper_ptr->label(), right_symbol,
                    parent_oper_ptr->label(), fma_other_arg_symbol);
                nflops_total += 2;  // extra flop due to accumulation + extra
                                    // flop due to FMA
              } else {
                os << context->accumulate_binary_expr(
                    curr_symbol, left_symbol, oper_ptr->label(), right_symbol);
                nflops_total += 1;  // extra flop due to accumulation
              }

            } else {  // assign, not accumulate
              if (generate_fma) {
                os << context->assign_ternary_expr(
                    parent_symbol, left_symbol, oper_ptr->label(), right_symbol,
                    parent_oper_ptr->label(), fma_other_arg_symbol);
                nflops_total += 1;  // extra flop due to FMA
              } else {
                os << context->assign_binary_expr(
                    curr_symbol, left_symbol, oper_ptr->label(), right_symbol);
              }
            }

            nflops_total += (1 + nflops(left_symbol) + nflops(right_symbol));
          }
          if (vectorize_by_line) os << line_vloop->close();

          // if produced FMA, do not forget to mark the parent scheduled
          if (generate_fma) {
            parent_oper_ptr->schedule();
          }

          goto next;
        }
      }

      // print simple assignment statement
      if (current_vertex->num_exit_arcs() == 1) {
        typedef DGArcDirect arc_type;
        std::shared_ptr<arc_type> arc_ptr =
            std::dynamic_pointer_cast<arc_type, DGArc>(
                *(current_vertex->first_exit_arc()));
        if (arc_ptr) {
#if CHECK_SAFETY
          current_vertex->declared(true);

          std::shared_ptr<DGVertex> rhs_arg = arc_ptr->dest();
          if (!rhs_arg->declared() && rhs_arg->need_to_compute()) {
            std::cout << "Current vertex:" << endl;
            current_vertex->print(std::cout);
            std::cout << "rhs_arg      :" << endl;
            rhs_arg->print(std::cout);
            throw ProgrammingError(
                "DirectedGraph::print_def() -- rhs_arg not declared");
          }
#endif

          // If this is an Integral in a target IntegralSet AND
          // can accumulate targets directly -- use '+=' instead of '='
          const bool accumulate_not_assign =
              accumulate_targets_directly &&
              IntegralInTargetIntegralSet()(current_vertex);

          if (context->comments_on()) {
            oss.str(null_str);
            oss << current_vertex->label()
                << (accumulate_not_assign ? " += " : " = ")
                << arc_ptr->dest()->label();
            os << context->comment(oss.str()) << endl;
          }

          // convert symbols to their vector form if needed
          std::string curr_symbol = current_vertex->symbol();
          std::string rhs_symbol = arc_ptr->dest()->symbol();
          if (vectorize) {
            curr_symbol = to_vector_symbol(current_vertex);
            rhs_symbol = to_vector_symbol(arc_ptr->dest());
          }

          if (vectorize_by_line) os << line_vloop->open();
          if (accumulate_not_assign) {
            os << context->accumulate(curr_symbol, rhs_symbol);
            nflops_total += nflops(rhs_symbol) + 1;  // +1 due to +=
          } else {
            os << context->assign(curr_symbol, rhs_symbol);
            nflops_total += nflops(rhs_symbol);
          }
          if (vectorize_by_line) os << line_vloop->close();

          goto next;
        }
      }

      // print out a recurrence relation
      if (current_vertex->num_exit_arcs() != 0) {
        // printing a recurrence relation
        // std::cout << "DirectedGraph::print_def(): a RR making " <<
        // current_vertex->description() << std::endl;
        typedef DGArcRR arc_type;
        std::shared_ptr<arc_type> arc_ptr =
            std::dynamic_pointer_cast<arc_type, DGArc>(
                *(current_vertex->first_exit_arc()));
        if (arc_ptr) {
          std::shared_ptr<RecurrenceRelation> rr = arc_ptr->rr();
          os << rr->spfunction_call(context, dims);
          nflops_total += rr->nflops();

          goto next;
        }
      }

#if 0
      {
        current_vertex->print(std::cout);
        throw std::runtime_error(
                                 "DirectedGraph::print_def() -- cannot handle this vertex yet");
      }
#endif
    }

  next:
    current_vertex->schedule();
    current_vertex = current_vertex->postcalc();

  } while (current_vertex != 0);

  os << outer_vloop->close();
  os << lsi_loop->close();
  os << hsi_loop->close();

  //
  // Accumulate targets
  //
  if (accumulate_targets_indirectly) {
    os << context->comment("Accumulate target integral sets") << std::endl;
    // const std::string& stack_name = registry()->stack_name();
    const bool vecdim_is_static = dims->vecdim_is_static();
#if 0
    unsigned int vecdim_rank;
    if (vecdim_is_static) {
      std::shared_ptr<CTimeEntity<int> > cptr = std::dynamic_pointer_cast<CTimeEntity<int> ,Entity >(dims->vecdim());
      vecdim_rank = cptr->value();
    }
    const std::string times_vecdim("*" + dims->vecdim_label());
#endif

    // Loop over targets
    unsigned int curr_target = 0;
    for (target_iter t = targets_.begin(); t != targets_.end();
         ++t, ++curr_target) {
      const ver_ptr& tptr = vertex_ptr(*t);

      size s = (tptr)->size();
      std::shared_ptr<CTimeEntity<int> > bdim(new CTimeEntity<int>(s));

      std::shared_ptr<Entity> bvecdim;
      if (!vecdim_is_static) {
        std::shared_ptr<RTimeEntity<EntityTypes::Int> > vecdim =
            std::dynamic_pointer_cast<RTimeEntity<EntityTypes::Int>, Entity>(
                dims->vecdim());
        bvecdim = vecdim * bdim;
      } else {
        std::shared_ptr<CTimeEntity<int> > vecdim =
            std::dynamic_pointer_cast<CTimeEntity<int>, Entity>(dims->vecdim());
        bvecdim = vecdim * bdim;
      }

      // For now write an explicit loop for each target. In the future should:
      // 1) check if all computed targets are adjacent (accumulated targets are
      // adjacent by the virtue of the allocation mechanism) 2) check the sizes
      // and insert optimized calls, if possible

      // form a single loop over integrals and vector dimension
      // NOTE: single loop suffices because if outer/inner strides are not 1
      // this block of code should not be executed

#if 0
      std::shared_ptr<Entity> loopmax;
      if (vecdim_is_static) {
const int loopmax_value = s*vecdim_rank;
loopmax = std::shared_ptr<Entity>(new CTimeEntity<int>(loopmax_value));
      }
      else {
ostringstream oss;  oss << s << times_vecdim;
loopmax = std::shared_ptr<Entity>(new RTimeEntity<EntityTypes::Int>(oss.str()));
      }

      std::string loopvar("i");
      ForLoop loop(context,loopvar,loopmax,ctimeconst_zero);
      os << loop.open();
      std::string acctarget;
      {
ostringstream oss;
oss << stack_name << "[" << target_accums_[curr_target] << times_vecdim << "+" << loopvar << "]";
acctarget = oss.str();
      }
      std::string target;
      {
ostringstream oss;
oss << stack_name << "[" << (tptr)->address() << times_vecdim << "+" << loopvar << "]";
target = oss.str();
      }
      os << context->accumulate(acctarget,target);
      os << loop.close();
#endif
      os << "_libint2_static_api_inc1_short_(" << registry()->stack_name()
         << "+" << target_accums_[curr_target] << "*" << dims->vecdim()->id()
         << "," << registry()->stack_name() << "+" << (tptr)->address() << "*"
         << dims->vecdim()->id() << "," << bvecdim->id() << ")" << endl;

      nflops_total += s;
    }
  }

  // Outside of loops stack symbols don't make sense, so we must define loop
  // variables hsi, lsi, and vi to 0
  os << context->decldef(context->type_name<const int>(), "hsi", "0");
  os << context->decldef(context->type_name<const int>(), "lsi", "0");
  os << context->decldef(context->type_name<const int>(), "vi", "0");

  //
  // Now pass back all targets through the inteval object, if needed.
  //
  if (registry()->return_targets()) {
    unsigned int curr_target = 0;
    for (target_iter t = targets_.begin(); t != targets_.end();
         ++t, ++curr_target) {
      const ver_ptr& tptr = vertex_ptr(*t);
      const std::string& symbol =
          (accumulate_targets_indirectly
               //                                                                    is this correct?         ???
               ? stack_symbol(context, target_accums_[curr_target],
                              (tptr)->size(), dims->low_label(),
                              dims->vecdim_label(), registry()->stack_name())
               : (tptr)->symbol());
      os << "inteval->targets[" << curr_target
         << "] = " << context->value_to_pointer(symbol)
         << context->end_of_stat() << endl;
    }
  }

  // Print out the number of flops
  oss.str(null_str);
  oss << "Number of flops = " << nflops_total;
  os << context->comment(oss.str()) << endl;

  if (context->cparams()->count_flops()) {
    oss.str(null_str);
    oss << nflops_total << " * " << dims->high_label() << " * "
        << dims->low_label() << " * " << dims->vecdim_label();
    os << context->assign_binary_expr("inteval->nflops[0]",
                                      "inteval->nflops[0]", "+", oss.str());
  }
}

void DirectedGraph::update_func_names() {
  // Loop over all vertices
  typedef vertices::iterator iter;
  for (iter v = stack_.begin(); v != stack_.end(); ++v) {
    const ver_ptr& vptr = vertex_ptr(*v);
    // for every vertex with children
    if ((vptr)->num_exit_arcs() > 0) {
      // if it must be computed using a RR
      std::shared_ptr<DGArc> arc = *((vptr)->first_exit_arc());
      std::shared_ptr<DGArcRR> arcrr =
          std::dynamic_pointer_cast<DGArcRR, DGArc>(arc);
      if (arcrr != 0) {
        std::shared_ptr<RecurrenceRelation> rr = arcrr->rr();
        // and the RR is complex (i.e. likely to result in a function call)
        if (!rr->is_simple()) {
          // add it to the RRStack
          func_names_[rr->label()] = true;
        }
      }
    }
  }
}

bool DirectedGraph::cannot_enclose_in_outer_vloop() const {
  std::shared_ptr<DGVertex> current_vertex = first_to_compute_;
  do {
    const int nchildren = current_vertex->num_exit_arcs();
    if (nchildren > 0) {
      arc_ptr aptr = *(current_vertex->first_exit_arc());
      std::shared_ptr<DGArcRR> aptr_cast =
          std::dynamic_pointer_cast<DGArcRR, arc>(aptr);
      // if this is a RR
      if (aptr_cast != 0) {
        // and a non-trivial one
        std::shared_ptr<RecurrenceRelation> rr = aptr_cast->rr();
        if (!rr->is_simple()) {
          // if this is an Uncontract_Integral call, return false
          std::shared_ptr<Uncontract_Integral_base> rr_ucb_ptr =
              std::dynamic_pointer_cast<Uncontract_Integral_base,
                                        RecurrenceRelation>(rr);
          if (rr_ucb_ptr)
            return false;
          else
            return true;
        }
      }
    }
    current_vertex = current_vertex->postcalc();
  } while (current_vertex != 0);

  return false;
}

void DirectedGraph::find_subtrees() {
  // need to condense expressions?
  if (!registry()->condense_expr()) return;

  // Find subtrees by starting from the targets and moving down ...
  typedef vertices::iterator iter;
  for (iter v = stack_.begin(); v != stack_.end(); ++v) {
    const ver_ptr& vptr = vertex_ptr(*v);
    if ((vptr)->is_a_target() && (vptr)->num_entry_arcs() == 0) {
      find_subtrees_from(vptr);
    }
  }
}

void DirectedGraph::find_subtrees_from(const std::shared_ptr<DGVertex>& v) {
  // is not on a subtree already
  if (!v->subtree()) {
    bool useless_subtree = false;

    //
    // Subtrees are useless in the following cases:
    // 1) root is computed via RRs, not explicitly
    //
    {
      if (v->num_exit_arcs() > 0) {
        std::shared_ptr<DGArc> arc = *(v->first_exit_arc());
        std::shared_ptr<DGArcRR> arc_rr =
            std::dynamic_pointer_cast<DGArcRR, DGArc>(arc);
        if (arc_rr) useless_subtree = true;
      }
    }

    // create subtree
    if (!useless_subtree) {
      std::shared_ptr<DRTree> stree = DRTree::CreateRootedAt(v);

      // Remove all trivial subtrees
      if (stree) {
#if DISABLE_SUBTREES
        if (stree->nvertices() >= 0) {
          stree->detach();
        }
#else
        if (stree->nvertices() < 3) {
          stree->detach();
        }
#endif
      }
    }

    // move on to children
    typedef DGVertex::ArcSetType::const_iterator aciter;
    const aciter abegin = v->first_exit_arc();
    const aciter aend = v->plast_exit_arc();
    for (aciter a = abegin; a != aend; ++a) {
      find_subtrees_from((*a)->dest());
    }
  }
}

namespace {
struct PrerequisiteNotComputed {
  bool operator()(const DirectedGraph::vertices::value_type& v) {
    return !v.second->precomputed() && v.second->num_exit_arcs() == 0;
  }
};
}  // namespace

/// return true if there are vertices with 0 children but not pre-computed
bool DirectedGraph::missing_prerequisites() const {
  bool missing_prereqs = false;
  if (!this->registry()->ignore_missing_prereqs()) {
#if 0
  missing_prereqs =
      find_if(this->stack_.begin(), this->stack_.end(), [](const vertices::value_type& v) {
                return v.second->precomputed() == false && v.second->num_exit_arcs() == 0;
             }) != this->stack_.end();
#else
    PrerequisiteNotComputed pred;
    missing_prereqs =
        std::find_if(stack_.begin(), stack_.end(), pred) != stack_.end();
#endif
  }
  return missing_prereqs;
}

////

namespace libint2 {

namespace {
#if USE_ASSOCCONTAINER_BASED_DIRECTEDGRAPH
typedef DirectedGraph::targets::value_type value_type;
struct __NotUnrolledIntegralSet {
  bool operator()(const value_type& v) { return NotUnrolledIntegralSet()(v); }
};
#endif
};  // namespace

bool nonunrolled_targets(const DirectedGraph::targets& targets) {
  typedef DirectedGraph::target_citer citer;
  citer end = targets.end();
#if USE_ASSOCCONTAINER_BASED_DIRECTEDGRAPH
  if (end != find_if(targets.begin(), end, __NotUnrolledIntegralSet()))
#else
  if (end != find_if(targets.begin(), end, NotUnrolledIntegralSet()))
#endif
    return true;
  else
    return false;
}

void extract_symbols(const std::shared_ptr<DirectedGraph>& dg) {
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  // symbol extractor
  {
    std::shared_ptr<ExtractExternSymbols> extractor(new ExtractExternSymbols);
    dg->foreach (*extractor);
    const ExtractExternSymbols::Symbols& symbols = extractor->symbols();
    // pass on to the symbol maintainer of the current task
    taskmgr.current().symbols()->add(symbols);
#if DEBUG
    // print out the symbols
    std::cout << "Recovered symbols from DirectedGraph for " << dg << std::endl;
    typedef ExtractExternSymbols::Symbols::const_iterator citer;
    citer end = symbols.end();
    for (citer t = symbols.begin(); t != end; ++t) std::cout << *t << std::endl;
#endif
  }
  // RR extractor
  {
    std::shared_ptr<ExtractRR> extractor(new ExtractRR);
    dg->foreach (*extractor);
    const ExtractRR::RRList& rrlist = extractor->rrlist();
    // pass on to the symbol maintainer of the current task
    taskmgr.current().symbols()->add(rrlist);
  }
}

//////////////////////
void PrerequisitesExtractor::operator()(const std::shared_ptr<DGVertex>& v) {
#if DEBUG
  std::cout << "PrerequisitesExtractor: considering " << v->description()
            << std::endl;
  v->print(std::cout);
#endif
  if (!v->precomputed() && v->num_exit_arcs() == 0) {
#if DEBUG
    std::cout << "PrerequisitesExtractor: found candidate " << v->description()
              << std::endl;
#endif

#define EXTRACTINTEGRALSETS 1
#if EXTRACTINTEGRALSETS
    // if this is an integral that was a member of a shell set, add the parent
    // set to the prerequisite level
    bool member_of_shellset = false;
    std::shared_ptr<DGVertex> parent_shellset;
    typedef DGVertex::ArcSetType::const_iterator citer;
    for (citer i = v->entry_arcs().begin();
         i != v->entry_arcs().end() && !member_of_shellset; ++i) {
      const std::shared_ptr<DGArc>& arc = *i;
      std::shared_ptr<DGArcRR> arc_cast =
          std::dynamic_pointer_cast<DGArcRR, DGArc>(arc);
      if (arc_cast) {
        std::shared_ptr<RecurrenceRelation> rr = arc_cast->rr();
        std::shared_ptr<IntegralSet_to_Integrals_base> rr_cast =
            std::dynamic_pointer_cast<IntegralSet_to_Integrals_base,
                                      RecurrenceRelation>(rr);
        if (rr_cast) {
          member_of_shellset = true;
          parent_shellset = rr->rr_target();
        }
      }
    }
    if (member_of_shellset) {
#if DEBUG
      std::cout << "PrerequisitesExtractor: " << v->description()
                << " is a member of a shell set, will add that instead"
                << std::endl;
#endif
      if (vertices.end() ==
          find(vertices.begin(), vertices.end(), parent_shellset)) {
        vertices.push_back(parent_shellset);
#if DEBUG
        std::cout << "PrerequisitesExtractor: extracted "
                  << parent_shellset->description() << std::endl;
#endif
      } else {
#if DEBUG
        std::cout
            << "PrerequisitesExtractor: candidate's parent already extracted"
            << std::endl;
#endif
      }
    } else {
      vertices.push_back(v);
#if DEBUG
      std::cout << "PrerequisitesExtractor: extracted " << v->description()
                << std::endl;
#endif
    }

#else

    vertices.push_back(v);
#if DEBUG
    std::cout << "PrerequisitesExtractor: extracted " << v->description()
              << std::endl;
#endif

#endif
  }
}
//////////////////////
void VertexPrinter::operator()(const std::shared_ptr<DGVertex>& v) {
  os << "VertexPrinter: " << v->description() << std::endl;
#if DEBUG
  v->print(os);
#endif
}

};  // namespace libint2
