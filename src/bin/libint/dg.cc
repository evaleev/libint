
#define USE_CG1 1

#include <rr.h>
#include <dg.h>
#include <strategy.h>

using namespace std;
using namespace libint2;

DGArc::DGArc(const SafePtr<DGVertex>& orig, const SafePtr<DGVertex>& dest) :
  orig_(orig), dest_(dest) {}

DGArcRR::DGArcRR(const SafePtr<DGVertex>& orig, const SafePtr<DGVertex>& dest) :
  DGArc(orig,dest) {}

DGVertex::DGVertex() :
  parents_(), children_(), target_(false), can_add_arcs_(true), num_tagged_arcs_(0),
  precalc_(), postcalc_()
{
}

DGVertex::DGVertex(const vector< SafePtr<DGArc> >& parents, const vector< SafePtr<DGArc> >& children) :
  parents_(parents), children_(children), target_(false), can_add_arcs_(true),
  num_tagged_arcs_(0), precalc_(), postcalc_()
{
}

DGVertex::~DGVertex()
{
}

void
DGVertex::make_a_target()
{
  target_ = true;
}

void
DGVertex::add_exit_arc(const SafePtr<DGArc>& arc)
{
  if (can_add_arcs_) {
    children_.push_back(arc);
    arc->dest()->add_entry_arc(arc);
  }
  else
    throw CannotAddArc("DGVertex::add_entry_arc() -- cannot add arcs anymore");
}

void
DGVertex::del_exit_arc(const SafePtr<DGArc>& arc)
{
  typedef vector< SafePtr<DGArc> > vectype;

  if (can_add_arcs_) {
    vectype::iterator pos = find(children_.begin(),children_.end(), arc);
    if (pos == children_.end()) {
      arc->dest()->del_entry_arc(arc);
      children_.erase(pos);
    }
    else
      return;
  }
  else
    throw CannotAddArc("DGVertex::del_entry_arc() -- cannot add/remove arcs anymore");
}

void
DGVertex::add_entry_arc(const SafePtr<DGArc>& arc)
{
  if (can_add_arcs_)
    parents_.push_back(arc);
  else
    throw CannotAddArc("DGVertex::add_entry_arc() -- cannot add arcs anymore");
}

void
DGVertex::del_entry_arc(const SafePtr<DGArc>& arc)
{
  vector< SafePtr<DGArc> >::iterator location = find(parents_.begin(), parents_.end(), arc);
  parents_.erase(location);
}

void
DGVertex::prepare_to_traverse()
{
  can_add_arcs_ = false;
  num_tagged_arcs_ = 0;
}

const unsigned int
DGVertex::tag()
{
  return ++num_tagged_arcs_;
}

const unsigned int
DGVertex::num_entry_arcs() const
{
  return parents_.size();
}

const unsigned int
DGVertex::num_exit_arcs() const
{
  return children_.size();
}

SafePtr<DGArc>
DGVertex::entry_arc(unsigned int p) const
{
  return parents_.at(p);
}

SafePtr<DGArc>
DGVertex::exit_arc(unsigned int c) const
{
  return children_.at(c);
}

void
DGVertex::reset()
{
  unsigned int nchildren = children_.size();
  for(int c=0; c<nchildren; c++) {
    children_[c]->dest()->del_entry_arc(children_[c]);
    children_[c].reset();
  }
  children_.resize(0);
  target_ = false;
  can_add_arcs_ = true;
  num_tagged_arcs_ = 0;
  precalc_.reset();
  postcalc_.reset();
}

void
DGVertex::set_label(const std::string& label)
{
  label_ = label;
}

///////////////////////////////////////////////////

DirectedGraph::DirectedGraph() :
  stack_(default_size_,SafePtr<DGVertex>()), first_free_(0), first_to_compute_()
{
}

DirectedGraph::~DirectedGraph()
{
}

void
DirectedGraph::append_target(const SafePtr<DGVertex>& target)
{
  target->make_a_target();
  try {
    add_vertex(target);
  }
  catch (VertexAlreadyOnStack) {
    return;
  }
}

bool
DirectedGraph::add_vertex(const SafePtr<DGVertex>& vertex)
{
  bool already_on_stack = false;
  for(int i=0; i<first_free_; i++) {
    if(vertex->equiv(stack_[i])) {
      already_on_stack = true;
      //break;
      return false;
    }
  }
  if(!already_on_stack) {
    if (first_free_ == stack_.size()) {
      stack_.resize( stack_.size() + default_size_ );
      cout << "Increased size of DirectedGraph's stack to "
           << stack_.size() << endl;
    }
    char label[80];  sprintf(label,"vertex%d",first_free_);
    vertex->set_label(label);
    stack_[first_free_++] = vertex;
    return true;
  }
  else
    throw VertexAlreadyOnStack("DirectedGraph::add_vertex() -- vertex already on stack");

}

bool
DirectedGraph::vertex_is_on(const SafePtr<DGVertex>& vertex) const
{
  for(int i=0; i<first_free_; i++)
    if(vertex->equiv(stack_[i]))
      return true;

  return false;
}

void
DirectedGraph::prepare_to_traverse()
{
  for(int i=0; i<first_free_; i++)
    stack_[i]->prepare_to_traverse();
}

void
DirectedGraph::traverse()
{
  // Initialization
  prepare_to_traverse();

  // Start at the targets which don't have parents
  for(int i=0; i<first_free_; i++) {
    if (stack_[i]->is_a_target() && stack_[i]->num_entry_arcs() == 0) {
      SafePtr<DGVertex> vertex_ptr = stack_[i];
      // First, since this target doesn't have parents we can schedule its computation
      schedule_computation(vertex_ptr);
      int nchildren = vertex_ptr->num_exit_arcs();
      for(int c=0; c<nchildren; c++)
        traverse_from(vertex_ptr->exit_arc(c));
    }
  }
}

void
DirectedGraph::traverse_from(const SafePtr<DGArc>& arc)
{
  SafePtr<DGVertex> orig = arc->orig();
  SafePtr<DGVertex> dest = arc->dest();
  if (dest->precomputed())
    return;
  const unsigned int num_tags = dest->tag();
  const unsigned int num_parents = dest->num_entry_arcs();
  /*if (num_tags == 1) {
    int nchildren = dest->num_exit_arcs();
    for(int c=0; c<nchildren; c++)
      traverse_from(dest->exit_arc(c));
      }*/
  if (num_tags == num_parents) {
    schedule_computation(dest);
    int nchildren = dest->num_exit_arcs();
    for(int c=0; c<nchildren; c++)
      traverse_from(dest->exit_arc(c));
  }
}

void
DirectedGraph::schedule_computation(const SafePtr<DGVertex>& vertex)
{
  vertex->set_precalc(SafePtr<DGVertex>());
  vertex->set_postcalc(first_to_compute_);
  if (first_to_compute_ != 0)
    first_to_compute_->set_precalc(vertex);
  first_to_compute_ = vertex;
}


void
DirectedGraph::debug_print_traversal(std::ostream& os) const
{
  SafePtr<DGVertex> current_vertex = first_to_compute_;

  os << "Debug print of traversal order" << endl;

  do {
    current_vertex->print(os);
    os << endl;
    current_vertex = current_vertex->postcalc();
  } while (current_vertex != 0);
}

void
DirectedGraph::print_to_dot(std::ostream& os) const
{
  os << "digraph G {" << endl
     << "  size = \"8,8\"" << endl;
  for(int i=0; i<first_free_; i++) {
    SafePtr<DGVertex> vertex = stack_[i];
    os << "  " << vertex->label()
       << " [ label = \"";
    vertex->print(os);
    os << "\"]" << endl;
  }
  for(int i=0; i<first_free_; i++) {
    SafePtr<DGVertex> vertex = stack_[i];
    unsigned int narcs = vertex->num_exit_arcs();
    for(int a=0; a<narcs; a++) {
      SafePtr<DGVertex> dest = vertex->exit_arc(a)->dest();
      os << "  " << vertex->label() << " -> "
         << dest->label() << endl;
    }
  }
  os << "}" << endl;
}

void
DirectedGraph::reset()
{
  // Reset each vertex, releasing all arcs
  for(int i=0; i<first_free_; i++)
    stack_[i]->reset();
  for(int i=0; i<first_free_; i++)
    stack_[i].reset();
  // if everything went OK then resize stack_ to 0
  stack_.resize(default_size_);
  first_free_ = 0;
  first_to_compute_.reset();
}



#if USE_GC1

/// Apply a strategy to all vertices not yet computed (i.e. which do not have exit arcs)
void
DirectedGraph::apply(const SafePtr<Strategy>& strategy)
{
  const int num_vertices_on_graph = first_free_;
  for(int v=0; v<num_vertices_on_graph; v++) {
    if (stack_[v]->num_exit_arcs() != 0)
      continue;

    SafePtr<DirectedGraph> this_ptr = SafePtr_from_this();
    SafePtr<RecurrenceRelation> rr0 = strategy->optimal_rr(this_ptr,stack_[v]);
    if (rr0 == 0)
      return;

    // add children to the graph
    SafePtr<DGVertex> target = rr0->rr_target();
    const int num_children = rr0->num_children();
    for(int c=0; c<num_children; c++) {
      SafePtr<DGVertex> child = rr0->rr_child(c);
      SafePtr<DGArc> arc(new DGArcRel<RecurrenceRelation>(target,child,rr0));
      target->add_exit_arc(arc);
    }

    // and apply strategy to the children
    for(int c=0; c<num_children; c++) {
      SafePtr<DGVertex> child = rr0->rr_child(c);
      apply_to(child,strategy);
    }
  }
}

/// Add vertex to graph and apply a strategy to vertex recursively
void
DirectedGraph::apply_to(const SafePtr<DGVertex>& vertex, const SafePtr<Strategy>& strategy)
{
  try {
    add_vertex(vertex);
  }
  catch (VertexAlreadyOnStack) {
    return;
  }

  SafePtr<RecurrenceRelation> rr0 = strategy->optimal_rr(SafePtr_from_this(),vertex);
  if (rr0 == 0)
    return;

  SafePtr<DGVertex> target = rr0->rr_target();
  const int num_children = rr0->num_children();
  for(int c=0; c<num_children; c++) {
    SafePtr<DGVertex> child = rr0->rr_child(c);
    SafePtr<DGArc> arc(new DGArcRel<RecurrenceRelation>(target,child,rr0));
    target->add_exit_arc(arc);
  }

  for(int c=0; c<num_children; c++) {
    SafePtr<DGVertex> child = rr0->rr_child(c);
    apply_to(child,strategy);
  }
}

#else

/// Apply a strategy to all vertices not yet computed (i.e. which do not have exit arcs)
void
DirectedGraph::apply(const SafePtr<Strategy>& strategy)
{
  const int num_vertices_on_graph = first_free_;
  for(int v=0; v<num_vertices_on_graph; v++) {
    if (stack_[v]->num_exit_arcs() != 0)
      continue;

    SafePtr<DirectedGraph> this_ptr = SafePtr_from_this();
    SafePtr<RecurrenceRelation> rr0 = strategy->optimal_rr(this_ptr,stack_[v]);
    if (rr0 == 0)
      return;

    // add children to the graph
    SafePtr<DGVertex> target = rr0->rr_target();
    const int num_children = rr0->num_children();
    for(int c=0; c<num_children; c++) {
      SafePtr<DGVertex> child = rr0->rr_child(c);
      SafePtr<DGArc> arc(new DGArcRel<RecurrenceRelation>(target,child,rr0));
      target->add_exit_arc(arc);
      apply_to(child,strategy);
    }

  }
}

/// Add vertex to graph and apply a strategy to vertex recursively
void
DirectedGraph::apply_to(const SafePtr<DGVertex>& vertex, const SafePtr<Strategy>& strategy)
{
  try {
    add_vertex(vertex);
  }
  catch (VertexAlreadyOnStack) {
    return;
  }

  SafePtr<RecurrenceRelation> rr0 = strategy->optimal_rr(SafePtr_from_this(),vertex);
  if (rr0 == 0)
    return;

  SafePtr<DGVertex> target = rr0->rr_target();
  const int num_children = rr0->num_children();
  for(int c=0; c<num_children; c++) {
    SafePtr<DGVertex> child = rr0->rr_child(c);
    SafePtr<DGArc> arc(new DGArcRel<RecurrenceRelation>(target,child,rr0));
    target->add_exit_arc(arc);
    apply_to(child,strategy);
  }

}


// Optimize out simple recurrence relations
void
DirectedGraph::optimize_rr_out()
{
  for(int v=0; v<first_free_; v++) {

    SafePtr<DGVertex> vertex = stack_[v];
    if (vertex->num_exit_arcs()) {
      SafePtr<DGArc> arc0 = vertex->exit_arc(0);
      SafePtr<DGArcRR> arc0_cast = dynamic_pointer_cast<DGArcRR,DGArc>(arc0);
      if (arc0_cast == 0)
        continue;
      SafePtr<RecurrenceRelation> rr = arc0_cast->rr();

      // Optimize if the recurrence relation is simple
      if (rr->is_simple()) {

        unsigned int nchildren = rr->num_children();
        unsigned int nexpr = rr->num_expr();

        cout << "RR: nchildren = " << nchildren << " nexpr = " << nexpr << endl;

        // Remove arcs connecting this vertex to children
        for(int c=0; c<nchildren; c++)
          vertex->del_exit_arc(vertex->exit_arc(c));

        // and instead insert expressions
        for(int e=0; e<nexpr; e++) {
          SafePtr< AlgebraicOperator<DGVertex> > rr_expr_cast = dynamic_pointer_cast<AlgebraicOperator<DGVertex>,DGVertex>(rr->rr_expr(e));
          if (rr_expr_cast)
            insert_expr_at(vertex,rr_expr_cast);
          else
            throw runtime_error("DirectedGraph::optimize_rr_out() -- expression of invalid type");
        }

      }
    }
  }
}


void
DirectedGraph::insert_expr_at(const SafePtr<DGVertex>& where, const SafePtr< AlgebraicOperator<DGVertex> >& expr)
{
  typedef AlgebraicOperator<DGVertex> ExprType;

  SafePtr<DGVertex> expr_vertex = dynamic_pointer_cast<DGVertex,ExprType>(expr);
  if (!add_vertex(expr_vertex))
    return;
  SafePtr<DGArc> arc(new DGArcDirect(where,expr_vertex));
  where->add_exit_arc(arc);

  // See if left operand is also an operator
  SafePtr<ExprType> left_cast = dynamic_pointer_cast<ExprType,DGVertex>(expr->left());
  if (left_cast)
    insert_expr_at(expr_vertex,left_cast);
  else {
    if (!add_vertex(expr->left()))
      return;
    SafePtr<DGArc> arc(new DGArcDirect(expr_vertex,expr->left()));
    expr_vertex->add_exit_arc(arc);
  }

  // See if right operand is also an operator
  SafePtr<ExprType> right_cast = dynamic_pointer_cast<ExprType,DGVertex>(expr->right());
  if (right_cast)
    insert_expr_at(expr_vertex,right_cast);
  else {
    if (!add_vertex(expr->right()))
      return;
    SafePtr<DGArc> arc(new DGArcDirect(expr_vertex,expr->right()));
    expr_vertex->add_exit_arc(arc);
  }
}

#endif

