
#include <algorithm>
#include <dgvertex.h>
#include <global_macros.h>

using namespace std;
using namespace libint2;

#define LOCAL_DEBUG 0
#define SAFE_DGVERTEX 0

DGVertex::DGVertex(ClassID tid) :
  dg_(0), subtree_(SafePtr<DRTree>()), typeid_(tid), parents_(), children_(), target_(false), can_add_arcs_(true), num_tagged_arcs_(0),
  postcalc_(), graph_label_(), referred_vertex_(0), nrefs_(0),
  symbol_(), address_(MemoryManager::InvalidAddress), need_to_compute_(true), instid_()
{
}

DGVertex::DGVertex(const DGVertex& v) :
  dg_(v.dg_), subtree_(v.subtree_), typeid_(v.typeid_), parents_(v.parents_), children_(v.children_), target_(v.target_),
  can_add_arcs_(v.can_add_arcs_), num_tagged_arcs_(v.num_tagged_arcs_),
  postcalc_(v.postcalc_), graph_label_(v.graph_label_),
  referred_vertex_(v.referred_vertex_), nrefs_(v.nrefs_),
  symbol_(v.symbol_), address_(v.address_), need_to_compute_(v.need_to_compute_), instid_(v.instid_)
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
#if SAFE_DGVERTEX
    typedef ArcSetType::const_iterator aciter;
    if (!children_.empty()) {
      const aciter abegin = children_.begin();
      const aciter aend = children_.end();
      for(aciter a=abegin; a!=aend; ++a) { 
	if ((*a)->dest() == child)
	  throw ProgrammingError("DGVertex::add_exit_arc() -- arc already exists");
      }
    }
#endif
    children_.push_back(arc);
    child->add_entry_arc(arc);
#if DEBUG
    std::cout << "add_exit_arc: added arc from " << arc->orig()->description() << " to " << arc->dest()->description() << std::endl;
#endif
  }
  else
    throw CannotAddArc("DGVertex::add_entry_arc() -- cannot add arcs anymore");
}

void
DGVertex::del_exit_arc(const SafePtr<DGArc>& arc)
{
  if (can_add_arcs_) {
    if (!children_.empty()) {
      ArcSetType::iterator pos = find(children_.begin(),children_.end(), arc);
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
  if (can_add_arcs_) {
    if (num_exit_arcs()) {
      do {
        del_exit_arc(*(children_.begin()));
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
      const aiter begin = children_.begin();
      const aiter end = children_.end();
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
    ArcSetType::iterator location = find(parents_.begin(), parents_.end(), arc);
    if (location != parents_.end())
      parents_.erase(location);
    else
      throw std::runtime_error("DGVertex::del_entry_arc() -- the arc doesn't exist");
  }
  else
    throw std::runtime_error("DGVertex::del_entry_arc() -- no arcs to delete");
}

void
DGVertex::detach()
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

namespace {
  struct __ArcDestEqual {
    __ArcDestEqual(const SafePtr<DGVertex>& v) : v_(v) {}
    bool operator()(const SafePtr<DGArc>& a) { return a->dest() == v_; }
    const SafePtr<DGVertex>& v_;
  };
}

const SafePtr<DGArc>&
DGVertex::exit_arc(const SafePtr<DGVertex>& v) const
{
  static SafePtr<DGArc> nullptr;
  __ArcDestEqual predicate(v);
  const ArcSetType::const_iterator end = children_.end();
  const ArcSetType::const_iterator pos = find_if(children_.begin(),children_.end(),predicate);
  if (pos != end)
    return *pos;
  else 
    return nullptr;
}

void
DGVertex::reset()
{
  dg_ = 0;
  subtree_ = SafePtr<DRTree>();

  typedef ArcSetType::const_iterator citer;
  typedef ArcSetType::iterator iter;
  const citer end = children_.end();
  for(iter a=children_.begin(); a!=end; ++a) {
    (*a)->dest()->del_entry_arc(*a);
    (*a).reset();
  }
  children_.clear();

  target_ = false;
  can_add_arcs_ = true;
  num_tagged_arcs_ = 0;
  postcalc_.reset();
  graph_label_.clear();
  reset_symbol();
  address_ = MemoryManager::InvalidAddress;
  need_to_compute_ = true;
  referred_vertex_ = 0;
  nrefs_ = 0;
}

const std::string&
DGVertex::graph_label() const
{
  if (!graph_label_.empty())
    return graph_label_;
  else
    throw GraphLabelNotSet("DGVertex::graph_label() -- graph label not set");
}

void
DGVertex::set_graph_label(const std::string& label)
{
  graph_label_ = label;
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
  referred_vertex_ = V.get();
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
DGVertex::symbol() const
{
  if (referred_vertex_ && referred_vertex_->symbol_set())
    return referred_vertex_->symbol();
  else {
    if (!symbol_.empty())
      return symbol_;
    else {
#if DEBUG
      cout << "DGVertex::symbol() -- symbol not set for " << description() << endl;
#endif
      throw SymbolNotSet("DGVertex::symbol() -- symbol not set");
    }
  }
}

void
DGVertex::set_symbol(const std::string& symbol)
{
  symbol_ = symbol;
#if DEBUG
  cout << "Set symbol for " << description() << " to " << symbol << endl;
#endif
}

void
DGVertex::reset_symbol()
{
  symbol_.clear();
}

DGVertex::Address
DGVertex::address() const
{
  if (referred_vertex_ && referred_vertex_->address_set())
    return referred_vertex_->address();
  else {
    if (address_ != MemoryManager::InvalidAddress)
      return address_;
    else {
      throw AddressNotSet("DGVertex::address() -- address not set");
    }
  }
}

void
DGVertex::set_address(const Address& address)
{
  address_ = address;
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

void
DGVertex::print(std::ostream& os) const
{
  using std::endl;
  std::string prefix("DGVertex::print: ");
  os << prefix << "label = " << label() << endl;
#if DEBUG || LOCAL_DEBUG
  os << prefix << "this = " << this << endl;
  if (referred_vertex_ != 0) {
    os << prefix << "refers_to = " << referred_vertex_ << endl;
  }
  else {
    if (symbol_set())
      os << prefix << "symbol = " << symbol() << endl;
    if (address_set())
      os << prefix << "address = " << address() << endl;
    os << prefix << "size = " << size() << endl;
  }
#endif
}

void
DGVertex::unregister() const
{
}

////

bool
UnrolledIntegralSet::operator()(const SafePtr<DGVertex>& V)
{
  const unsigned int outdegree = V->num_exit_arcs();
  if (outdegree == 0) return false;
  
  const SafePtr<DGArc> arc0 = *(V->first_exit_arc());
  // Is this DGArcRR?
  const SafePtr<DGArcRR> arcrr = dynamic_pointer_cast<DGArcRR,DGArc>(arc0);
  if (arcrr == 0) return false;
  // Is this DGArcRR<IntegralSet_to_Integral>? If invariant_type() is false, then yes
  return !arcrr->rr()->invariant_type();
}

bool
NotUnrolledIntegralSet::operator()(const SafePtr<DGVertex>& V)
{
  return !UnrolledIntegralSet()(V);
}

bool
IntegralInTargetIntegralSet::operator()(const SafePtr<DGVertex>& V)
{
  const unsigned int indegree = V->num_entry_arcs();
  if (indegree != 1) return false;
  const SafePtr<DGVertex>& parent = (*(V->first_entry_arc()))->orig();
  if (parent->is_a_target()) return UnrolledIntegralSet()(parent);
  return false;
}
