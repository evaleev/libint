
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
  parents_(0), children_(0), target_(false)
{
}

DGVertex::DGVertex(const vector<DGArc*>& parents, const vector<DGArc*>& children) :
  parents_(parents), children_(children), target_(false)
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
  children_.push_back(arc);
  arc->dest()->add_entry_arc(arc);
}

void
DGVertex::add_entry_arc(DGArc* arc)
{
  parents_.push_back(arc);
}


///////////////////////////////////////////////////

DirectedGraph::DirectedGraph() :
  stack_(default_size_), first_free_(0)
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


DGVertex*
DirectedGraph::traverse()
{
  // Initialization
  for(int i=0; i<first_free_; i++)
    stack_[i]->prepare_to_traverse();

  // Start at each target
  for(int i=0; i<first_free_; i++) {
    if (stack_[i]->is_a_target())
      DGVertex* vertex_ptr = stack_[i];
      int nchildren = vertex_ptr->num_children();
      for(int c=0; c<nchildren; c++)
        traverse_from(vertex_ptr,vertex_ptr->child(c));
  }

}
