
#include <algorithm>
#include <dgvertex.h>
#include <global_macros.h>

using namespace std;
using namespace libint2;

#define LOCAL_DEBUG 0

DGVertex::DGVertex(ClassID tid) :
  typeid_(tid), instid_(), dg_(0), graph_label_(), referred_vertex_(0),
  nrefs_(0), symbol_(), address_(MemoryManager::InvalidAddress), need_to_compute_(true),
#if CHECK_SAFETY
  declared_(false),
#endif
  parents_(), children_(), target_(false), can_add_arcs_(true), num_tagged_arcs_(0),
  postcalc_(), subtree_(SafePtr<DRTree>())
{
}

DGVertex::DGVertex(const DGVertex& v) :
  typeid_(v.typeid_), instid_(v.instid_), dg_(v.dg_), graph_label_(v.graph_label_), referred_vertex_(v.referred_vertex_),
  nrefs_(v.nrefs_), symbol_(v.symbol_), address_(v.address_), need_to_compute_(v.need_to_compute_),
#if CHECK_SAFETY
  declared_(v.declared_),
#endif
  parents_(v.parents_), children_(v.children_), target_(v.target_),
  can_add_arcs_(v.can_add_arcs_), num_tagged_arcs_(v.num_tagged_arcs_),
  postcalc_(v.postcalc_), subtree_(v.subtree_)
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

    // check if such arc exists already
    if (!children_.empty()) {
      typedef ArcSetType::const_iterator aciter;
      const aciter abegin = children_.begin();
      const aciter aend = children_.end();
      for(aciter a=abegin; a!=aend; ++a) {
	if ((*a)->dest() == child)
	  return;
      }
    }

    children_.push_back(arc);
    child->add_entry_arc(arc);
#if DEBUG
    std::cout << "add_exit_arc: added arc from " << arc->orig()->description() << " to " << arc->dest()->description() << std::endl;
#endif
  }
  else
    throw CannotAddArc("DGVertex::add_exit_arc() -- cannot add arcs anymore");
}

void
DGVertex::del_exit_arc(const SafePtr<DGArc>& arc)
{
  if (can_add_arcs_) {
    if (!children_.empty()) {
      ArcSetType::iterator pos = find(children_.begin(),children_.end(), arc);
      if (pos != children_.end()) {
        arc->dest()->del_entry_arc(arc);
#if DEBUG
    std::cout << "del_exit_arc: removed arc from " << arc->orig()->description() << " to " << arc->dest()->description() << std::endl;
#endif
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
#if CHECK_SAFETY
      aiter posB = find(begin,end,B);
      bool B_already_exists = (posB != end);
      if (B_already_exists)
        throw std::runtime_error("DGVertex::replace_exit_arc(A,B) -- arc B is found among children");
#endif
#if DEBUG || DEBUG_RESTRUCTURE
      std::cout << "replace_exit_arc: replacing arc from " << A->orig().get() << " to " << A->dest().get() << endl;
      std::cout << "replace_exit_arc:      with arc from " << B->orig().get() << " to " << B->dest().get() << endl;
#endif
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
  if (arc->orig() == arc->dest())
    throw CannotAddArc("DGVertex::add_entry_arc() -- arc connects node to itself");

  if (can_add_arcs_)
    parents_.push_back(arc);
  else
    throw CannotAddArc("DGVertex::add_entry_arc() -- cannot add arcs anymore");
#if DEBUG || DEBUG_RESTRUCTURE
  std::cout << "add_entry_arc: from " << arc->orig()->description() << " to " << arc->dest()->description() << std::endl;
  print(std::cout);
#endif
}

void
DGVertex::del_entry_arc(const SafePtr<DGArc>& arc)
{
  if (!parents_.empty()) {
    ArcSetType::iterator location = find(parents_.begin(), parents_.end(), arc);
    if (location != parents_.end()) {
      parents_.erase(location);
#if DEBUG || DEBUG_RESTRUCTURE
      std::cout << "del_entry_arc: removed arc from " << (*location)->orig()->description()
                << " to " << (*location)->dest()->description() << endl;
#endif
    }
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
  const unsigned int narcs = num_entry_arcs();
  if (narcs == 0)
    DGVertex::del_exit_arcs();
  else
    throw CannotPerformOperation("DGVertex::detach() -- cannot detach a vertex if it has entry arcs");
}

void
DGVertex::prepare_to_traverse()
{
  //can_add_arcs_ = false;
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
  static SafePtr<DGArc> nullptr_;
  __ArcDestEqual predicate(v);
  const ArcSetType::const_iterator end = children_.end();
  const ArcSetType::const_iterator pos = find_if(children_.begin(),children_.end(),predicate);
  if (pos != end)
    return *pos;
  else
    return nullptr_;
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
#if DEBUG
      cout << "DGVertex::refer_this_to() -- vertex " << description() << " will refer to " << V->description() << endl;
#endif
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
      if (referred_vertex_)
        cout << "DGVertex::symbol() -- referred_vertex_ = " << referred_vertex_->description() << endl;
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
  os << prefix << "this = " << this << endl;
  if (referred_vertex_ != 0) {
    os << prefix << "refers_to = " << referred_vertex_ << endl;
  }
  else {
    os << prefix << "precomputed = " << precomputed() << endl;
    if (symbol_set())
      os << prefix << "symbol = " << symbol() << endl;
    if (address_set())
      os << prefix << "address = " << address() << endl;
    os << prefix << "size = " << size() << endl;
    os << prefix << "next to compute = " << postcalc() << endl;
    os << prefix << "nparents = " << num_entry_arcs() << endl;
    unsigned int i=0;
    for(ArcSetType::const_iterator p=first_entry_arc(); p!=plast_entry_arc(); ++i,++p)
      os << prefix << "  parent " << i << ": " << (*p)->orig() << endl;
    os << prefix << "nchildren = " << num_exit_arcs() << endl;
    i=0;
    for(ArcSetType::const_iterator c=first_exit_arc(); c!=plast_exit_arc(); ++i,++c)
      os << prefix << "  child " << i << ": " << (*c)->dest() << endl;
    os << prefix << "ntags = " << num_tagged_arcs_ << endl;
  }
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
