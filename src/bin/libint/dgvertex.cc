
#include <dgvertex.h>

using namespace std;
using namespace libint2;

DGVertex::DGVertex() :
  parents_(), children_(), target_(false), can_add_arcs_(true), num_tagged_arcs_(0),
  precalc_(), postcalc_(), graph_label_(), referred_vertex_(SafePtr<DGVertex>()), nrefs_(0),
  symbol_(), address_(), need_to_compute_(true)
{
}

DGVertex::DGVertex(const vector< SafePtr<DGArc> >& parents, const vector< SafePtr<DGArc> >& children) :
  parents_(parents), children_(children), target_(false), can_add_arcs_(true),
  num_tagged_arcs_(0), precalc_(), postcalc_(), graph_label_(),
  referred_vertex_(SafePtr<DGVertex>()), nrefs_(0), symbol_(), address_(),
  need_to_compute_(true)
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
    SafePtr<DGVertex> child = arc->dest();
    const unsigned int nchildren = children_.size();
    for(int i = 0; i<nchildren; i++)
      if (children_[i]->dest() == child)
        return;
    children_.push_back(arc);
    child->add_entry_arc(arc);
  }
  else
    throw CannotAddArc("DGVertex::add_entry_arc() -- cannot add arcs anymore");
}

void
DGVertex::del_exit_arc(const SafePtr<DGArc>& arc)
{
  typedef vector< SafePtr<DGArc> > vectype;

  if (can_add_arcs_) {
    if (!children_.empty()) {
      vectype::iterator pos = find(children_.begin(),children_.end(), arc);
      if (pos != children_.end()) {
        arc->dest()->del_entry_arc(arc);
        children_.erase(pos);
      }
      else
        throw std::runtime_error("DGVertex::del_exit_arc() -- arc does not exist");
    }
    else
      throw std::runtime_error("DGVertex::del_exit_arc() -- no arcs to delete");
  }
  else
    throw CannotAddArc("DGVertex::del_exit_arc() -- cannot add/remove arcs anymore");
}

void
DGVertex::del_exit_arcs()
{
  typedef vector< SafePtr<DGArc> > vectype;

  if (can_add_arcs_) {
    if (num_exit_arcs()) {
      do {
        del_exit_arc(children_[0]);
      } while (num_exit_arcs() != 0);
    }
  }
  else
    throw CannotAddArc("DGVertex::del_exit_arcs() -- cannot add/remove arcs anymore");
}

void
DGVertex::replace_exit_arc(const SafePtr<DGArc>& A, const SafePtr<DGArc>& B)
{
  if (can_add_arcs_) {
    typedef ArcSetType::iterator aiter;
    if (!children_.empty()) {
      aiter begin = children_.begin();
      aiter end = children_.end();
      aiter posB = find(begin,end,B);
      bool B_already_exists = (posB != end);
      if (B_already_exists)
        throw std::runtime_error("DGVertex::replace_exit_arc(A,B) -- arc B is found among children");
      aiter posA = find(begin,end,A);
      if (posA != end) {
        *posA = B;
        A->dest()->del_entry_arc(A);
        B->dest()->add_entry_arc(B);
      }
      else
        throw std::runtime_error("DGVertex::replace_exit_arc(A,B) -- arc A is not found among exit arcs");
    }
    else
      throw CannotAddArc("DGVertex::replace_exit_arc() -- no arcs to replace");
  }
  else
    throw CannotAddArc("DGVertex::replace_exit_arc() -- cannot add/remove arcs anymore");
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
  if (!parents_.empty()) {
    vector< SafePtr<DGArc> >::iterator location = find(parents_.begin(), parents_.end(), arc);
    if (location != parents_.end())
      parents_.erase(location);
    else
      throw std::runtime_error("DGVertex::del_entry_arc() -- the arc doesn't exist");
  }
  else
    throw std::runtime_error("DGVertex::del_entry_arc() -- no arcs to delete");
}

void
DGVertex::detach() throw (CannotPerformOperation)
{
  // If there are no entry arcs -- then other vertices do not depend on this guy
  // Can safely remove exit arcs
  int narcs = num_entry_arcs();
  if (num_entry_arcs() == 0)
    DGVertex::del_exit_arcs();
  else
    throw CannotPerformOperation("DGVertex::detach() -- cannot detach a vertex if it has entry arcs");
}

void
DGVertex::prepare_to_traverse()
{
  can_add_arcs_ = false;
  num_tagged_arcs_ = 0;
}

unsigned int
DGVertex::tag()
{
  return ++num_tagged_arcs_;
}

unsigned int
DGVertex::num_entry_arcs() const
{
  return parents_.size();
}

unsigned int
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

SafePtr<DGArc>
DGVertex::exit_arc(const SafePtr<DGVertex>& v) const
{
  unsigned int nchildren = children_.size();
  for(int c=0; c<nchildren; c++) {
    SafePtr<DGArc> arc = children_[c];
    if (arc->dest() == v)
      return arc;
  }
  return SafePtr<DGArc>();
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
  graph_label_ = SafePtr<std::string>();
  reset_symbol();
  address_ = SafePtr<Address>();
  need_to_compute_ = true;
  referred_vertex_ = SafePtr<DGVertex>();
  nrefs_ = 0;
}

const std::string&
DGVertex::graph_label() const throw(GraphLabelNotSet)
{
  if (graph_label_)
    return *graph_label_;
  else
    throw GraphLabelNotSet("DGVertex::graph_label() -- graph label not set");
}

void
DGVertex::set_graph_label(const std::string& label)
{
  SafePtr<std::string> graph_label(new std::string(label));
  graph_label_ = graph_label;
}

void
DGVertex::refer_this_to(const SafePtr<DGVertex>& V)
{
  if (referred_vertex_ != 0) {
    if (referred_vertex_->equiv(V))
      return;
    else
      throw std::logic_error("DGVertex::refer_this_to() -- already referring to some other vertex");
  }
  cout << "Referring " << this->description() << " to " << V->description() << endl;
  referred_vertex_ = V;
  V->inc_nrefs();
}

void
DGVertex::inc_nrefs()
{
  if (nrefs_)
    throw std::logic_error("DGVertex::inc_nrefs() -- already referred to by another vertex");
  else
    ++nrefs_;
}

const std::string&
DGVertex::symbol() const throw(SymbolNotSet)
{
  if (referred_vertex_ && referred_vertex_->symbol_set())
    return referred_vertex_->symbol();
  else {
    if (symbol_)
      return *symbol_;
    else
      throw SymbolNotSet("DGVertex::symbol() -- symbol not set");
  }
}

void
DGVertex::set_symbol(const std::string& symbol)
{
  SafePtr<std::string> ptr(new std::string(symbol));
  symbol_ = ptr;
}

void
DGVertex::reset_symbol()
{
  symbol_ = SafePtr<std::string>();
}

DGVertex::Address
DGVertex::address() const throw(AddressNotSet)
{
  if (referred_vertex_ && referred_vertex_->address_set())
    return referred_vertex_->address();
  else {
    if (address_)
      return *address_;
    else {
      throw AddressNotSet("DGVertex::address() -- address not set");
    }
  }
}

void
DGVertex::set_address(const Address& address)
{
  SafePtr<Address> ptr(new Address(address));
  address_ = ptr;
}

void
DGVertex::need_to_compute(bool ntc)
{
  need_to_compute_ = ntc;
}

bool
DGVertex::need_to_compute() const
{
  if (referred_vertex_)
    return referred_vertex_->need_to_compute();
  else
    return need_to_compute_;
}

bool
DGVertex::precomputed() const
{
  if (referred_vertex_)
    return referred_vertex_->precomputed();
  else {
    return this_precomputed();
  }
}

