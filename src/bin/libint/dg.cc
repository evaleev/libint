
#include <rr.h>

using namespace std;
using namespace libint2;

DGArc::DGArc(const DGVertex* orig, const DGVertex* dest)
{
  orig_ = orig;
  dest_ = dest;
}

DGArc::~DGArc()
{
}

DGVertex::DGVertex() :
  parents_(0), children_(0)
{
}

DGVertex::DGVertex(const vector<DGArc*>& parents, const vector<DGArc*>& children) :
  parents_(parents), children_(children)
{
}

DGVertex::~DGVertex()
{
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
    stack_->~DGVertex;
}

/// Apply RR to target
template <class I, class RR>
void recurse<I,RR>(const I* target)
{
  add_vertex(target);
  
  RR* rr0 = new RR(target);
  const int num_children = rr0->num_children();

  for(int c=0; c<num_children; c++) {

    const DGVertex* child = rr0->children(c);
    add_vertex(child);

    DGArc* arc = DGArcRel<RR>(target,child,rr0);
    target->add_exit(arc);
    child->add_entry(arc);

    recurse<I,RR>(child);
    
  }

}
