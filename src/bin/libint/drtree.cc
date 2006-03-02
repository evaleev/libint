
#include <drtree.h>
#include <dgvertex.h>

#define LOCAL_DEBUG 0

using namespace libint2;

SafePtr<DRTree>
DRTree::CreateRootedAt(const SafePtr<DGVertex>& v)
{
  if (!v->subtree()) {
    SafePtr<DRTree> result(new DRTree(v));
    // Ugly that I have to add the vertex outside the constructor,
    // but enable_shared_from_this requires that a valid shared_ptr to this already exists
    result->grow();
  }
  else
    return v->subtree();
}

DRTree::DRTree(const SafePtr<DGVertex>& r) :
nvertices_(0), root_(r)
{
}

DRTree::~DRTree()
{
}

void
DRTree::grow()
{
  add_vertex(root());
}

const SafePtr<DGVertex>&
DRTree::root() const
{
  return root_;
}

void
DRTree::add_vertex(const SafePtr<DGVertex>& vertex)
{
  // If not root and has more than 1 parent -- it is not on the tree
  if ( (vertex->num_entry_arcs() <= 1 || vertex == root()) ) {
    if (vertex->subtree_)
      throw ProgrammingError("DRTree::add_vertex() -- vertex is on a subtree already");
    vertex->subtree_ = EnableSafePtrFromThis<this_type>::SafePtr_from_this();
    ++nvertices_;
#if LOCAL_DEBUG
    std::cout << "Vertex " << vertex->label() << " is on the following subtree:" << std::endl;
    std::cout << "  Root = " << root()->label() << std::endl;
    std::cout << "  nvertices = " << nvertices_ << std::endl;
#endif
    const unsigned int nchildren = vertex->num_exit_arcs();
    for(unsigned int c=0; c<nchildren; c++) {
      add_vertex(vertex->exit_arc(c)->dest());
    }
  }
}

void
DRTree::detach()
{
  detach_from(root());
}

void
DRTree::detach_from(const SafePtr<DGVertex>& v)
{
  if (v->subtree_.get() != this)
    return;
  else {
    v->subtree_ = SafePtr<DRTree>();
    const unsigned int nchildren = v->num_exit_arcs();
    for(unsigned int c=0; c<nchildren; c++)
      detach_from(v->exit_arc(c)->dest());
  }
}

