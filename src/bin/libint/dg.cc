
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
  int size = stack_.size();
  for(int i=0; i<size; i++)
    stack_[i]->~DGVertex();
}

void
DirectedGraph::add_vertex(DGVertex* vertex)
{
  stack_[first_free_++] = vertex;
}

