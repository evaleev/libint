
#include <rr.h>
#include <dg.h>

using namespace std;
using namespace libint2;

DGArc::DGArc(DGVertex* orig, DGVertex* dest)
{
  orig_ = orig;
  dest_ = dest;
}

DGArc::~DGArc()
{
}

DGVertex::DGVertex() :
  parents_(), children_(), target_(false), can_add_arcs_(true), num_tagged_arcs_(0),
  precalc_(0), postcalc_(0)
{
}

DGVertex::DGVertex(const vector<DGArc*>& parents, const vector<DGArc*>& children) :
  parents_(parents), children_(children), target_(false), can_add_arcs_(true),
  num_tagged_arcs_(0), precalc_(0), postcalc_(0)
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
DGVertex::add_exit_arc(DGArc* arc)
{
  if (can_add_arcs_) {
    children_.push_back(arc);
    arc->dest()->add_entry_arc(arc);
  }
  else
    throw CannotAddArc("DGVertex::add_entry_arc() -- cannot add arcs anymore");
}

void
DGVertex::add_entry_arc(DGArc* arc)
{
  if (can_add_arcs_)
    parents_.push_back(arc);
  else
    throw CannotAddArc("DGVertex::add_entry_arc() -- cannot add arcs anymore");
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

DGArc*
DGVertex::entry_arc(unsigned int p) const
{
  return parents_.at(p);
}

DGArc*
DGVertex::exit_arc(unsigned int c) const
{
  return children_.at(c);
}

///////////////////////////////////////////////////

DirectedGraph::DirectedGraph() :
  stack_(default_size_), first_free_(0), first_to_compute_(0)
{
}

DirectedGraph::~DirectedGraph()
{
  for(int i=0; i<first_free_; i++)
    if (!stack_[i]->is_a_target())
      stack_[i]->~DGVertex();
}

void
DirectedGraph::add_vertex(DGVertex* vertex)
{
  bool already_on_stack = false;
  for(int i=0; i<first_free_; i++) {
    if(vertex->equiv(stack_[i])) {
      already_on_stack = true;
      break;
    }
  }
  if(!already_on_stack) {
    if (first_free_ == stack_.size()) {
      stack_.resize( stack_.size() + default_size_ );
      cout << "Increased size of DirectedGraph's stack to "
           << stack_.size() << endl;
    }
    stack_[first_free_++] = vertex;
  }
  else
    throw VertexAlreadyOnStack("DirectedGraph::add_vertex() -- vertex already on stack");

}


void
DirectedGraph::prepare_to_traverse()
{
  for(int i=0; i<first_free_; i++)
    stack_[i]->prepare_to_traverse();
}

DGVertex*
DirectedGraph::traverse()
{
  // Initialization
  prepare_to_traverse();

  // Start at the targets which don't have parents
  for(int i=0; i<first_free_; i++) {
    if (stack_[i]->is_a_target() && stack_[i]->num_entry_arcs() == 0) {
      DGVertex* vertex_ptr = stack_[i];
      // First, since this target doesn't have parents we can schedule its computation
      schedule_computation(vertex_ptr);
      int nchildren = vertex_ptr->num_exit_arcs();
      for(int c=0; c<nchildren; c++)
        traverse_from(vertex_ptr->exit_arc(c));
    }
  }

}

void
DirectedGraph::traverse_from(DGArc* arc)
{
  DGVertex* orig = arc->orig();
  DGVertex* dest = arc->dest();
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
DirectedGraph::schedule_computation(DGVertex* vertex)
{
  vertex->set_precalc(0);
  vertex->set_postcalc(first_to_compute_);
  if (first_to_compute_ != 0)
    first_to_compute_->set_precalc(vertex);
  first_to_compute_ = vertex;
}


void
DirectedGraph::debug_print_traversal(std::ostream& os) const
{
  DGVertex* current_vertex = first_to_compute_;

  os << "Debug print of traversal order" << endl;

  do {
    current_vertex->print(os);
    current_vertex = current_vertex->postcalc();
  } while (current_vertex != 0);
}
