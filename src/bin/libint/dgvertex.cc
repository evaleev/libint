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

#include <dgvertex.h>
#include <global_macros.h>
#include <rr.h>

#include <algorithm>

using namespace std;
using namespace libint2;

#define LOCAL_DEBUG 0

DGVertex::DGVertex(ClassID tid)
    : typeid_(tid),
      instid_(),
      dg_(0),
      graph_label_(),
      referred_vertex_(0),
      refs_(),
      symbol_(),
      address_(MemoryManager::InvalidAddress),
      need_to_compute_(true),
#if CHECK_SAFETY
      declared_(false),
#endif
      parents_(),
      children_(),
      target_(false),
      can_add_arcs_(true),
      num_tagged_arcs_(0),
      postcalc_(),
      scheduled_(false),
      subtree_(std::shared_ptr<DRTree>()) {
}

DGVertex::DGVertex(const DGVertex& v)
    : typeid_(v.typeid_),
      instid_(v.instid_),
      dg_(v.dg_),
      graph_label_(v.graph_label_),
      referred_vertex_(v.referred_vertex_),
      refs_(v.refs_),
      symbol_(v.symbol_),
      address_(v.address_),
      need_to_compute_(v.need_to_compute_),
#if CHECK_SAFETY
      declared_(v.declared_),
#endif
      parents_(v.parents_),
      children_(v.children_),
      target_(v.target_),
      can_add_arcs_(v.can_add_arcs_),
      num_tagged_arcs_(v.num_tagged_arcs_),
      postcalc_(v.postcalc_),
      scheduled_(false),
      subtree_(v.subtree_) {
}

DGVertex::~DGVertex() {}

void DGVertex::make_a_target() { target_ = true; }

void DGVertex::add_exit_arc(const std::shared_ptr<DGArc>& arc) {
  if (can_add_arcs_) {
    std::shared_ptr<DGVertex> child = arc->dest();

    // check if such arc exists already
    if (!children_.empty()) {
      typedef ArcSetType::const_iterator aciter;
      const aciter abegin = children_.begin();
      const aciter aend = children_.end();
      for (aciter a = abegin; a != aend; ++a) {
        if ((*a)->dest() == child) return;
      }
    }

    children_.push_back(arc);
    child->add_entry_arc(arc);
#if DEBUG
    std::cout << "add_exit_arc: added arc " << arc << " from "
              << arc->orig()->description() << " to "
              << arc->dest()->description() << std::endl;
#endif
  } else
    throw CannotAddArc("DGVertex::add_exit_arc() -- cannot add arcs anymore");
}

void DGVertex::del_exit_arc(const std::shared_ptr<DGArc>& arc) {
  if (can_add_arcs_) {
    if (!children_.empty()) {
      ArcSetType::iterator pos = find(children_.begin(), children_.end(), arc);
      if (pos != children_.end()) {
        arc->dest()->del_entry_arc(arc);
#if DEBUG
        std::cout << "del_exit_arc: removed arc from "
                  << arc->orig()->description() << " to "
                  << arc->dest()->description() << std::endl;
#endif
        children_.erase(pos);
      } else
        throw std::runtime_error(
            "DGVertex::del_exit_arc() -- arc does not exist");
    } else
      throw std::runtime_error("DGVertex::del_exit_arc() -- no arcs to delete");
  } else
    throw CannotAddArc(
        "DGVertex::del_exit_arc() -- cannot add/remove arcs anymore");
}

void DGVertex::del_exit_arcs() {
  if (can_add_arcs_) {
    if (num_exit_arcs()) {
      do {
#if DEBUG_RESTRUCTURE
        std::cout << "DGVertex::del_exit_arcs(): num_exit_arcs = "
                  << this->num_exit_arcs() << std::endl;
#endif
#if DEBUG_RESTRUCTURE
        std::cout << "DGVertex::del_exit_arcs(): trying to delete exit arc: "
                  << children_.front().get() << std::endl;
        children_.front()->print(std::cout);
        std::cout.flush();
#endif
        del_exit_arc(*(children_.begin()));
#if DEBUG_RESTRUCTURE
        std::cout << "DGVertex::del_exit_arcs(): delete successful"
                  << std::endl;
#endif
      } while (num_exit_arcs() != 0);
    }
  } else
    throw CannotAddArc(
        "DGVertex::del_exit_arcs() -- cannot add/remove arcs anymore");
}

void DGVertex::replace_exit_arc(const std::shared_ptr<DGArc>& A,
                                const std::shared_ptr<DGArc>& B) {
  if (can_add_arcs_) {
    typedef ArcSetType::iterator aiter;
    if (!children_.empty()) {
      const aiter begin = children_.begin();
      const aiter end = children_.end();
#if CHECK_SAFETY
      aiter posB = find(begin, end, B);
      bool B_already_exists = (posB != end);
      if (B_already_exists)
        throw std::runtime_error(
            "DGVertex::replace_exit_arc(A,B) -- arc B is found among children");
#endif
#if DEBUG || DEBUG_RESTRUCTURE
      std::cout << "replace_exit_arc: replacing arc from " << A->orig().get()
                << " to " << A->dest().get() << endl;
      std::cout << "replace_exit_arc:      with arc from " << B->orig().get()
                << " to " << B->dest().get() << endl;
#endif
      aiter posA = find(begin, end, A);
      if (posA != end) {
        *posA = B;
        A->dest()->del_entry_arc(A);
        B->dest()->add_entry_arc(B);
      } else
        throw std::runtime_error(
            "DGVertex::replace_exit_arc(A,B) -- arc A is not found among exit "
            "arcs");
    } else
      throw CannotAddArc("DGVertex::replace_exit_arc() -- no arcs to replace");
  } else
    throw CannotAddArc(
        "DGVertex::replace_exit_arc() -- cannot add/remove arcs anymore");
}

void DGVertex::add_entry_arc(const std::shared_ptr<DGArc>& arc) {
  if (arc->orig() == arc->dest()) {
    std::cout << "DGVertex::add_entry_arc() : arc->orig = "
              << arc->orig()->description() << std::endl;
    std::cout << "DGVertex::add_entry_arc() : arc->dest = "
              << arc->dest()->description() << std::endl;
    throw CannotAddArc(
        "DGVertex::add_entry_arc() -- arc connects node to itself");
  }

  if (can_add_arcs_)
    parents_.push_back(arc);
  else
    throw CannotAddArc("DGVertex::add_entry_arc() -- cannot add arcs anymore");
#if DEBUG || DEBUG_RESTRUCTURE
  std::cout << "add_entry_arc: arc " << arc << " from "
            << arc->orig()->description() << " to "
            << arc->dest()->description() << std::endl;
  print(std::cout);
#endif
}

void DGVertex::del_entry_arc(const std::shared_ptr<DGArc>& arc) {
  if (!parents_.empty()) {
    ArcSetType::iterator location = find(parents_.begin(), parents_.end(), arc);
    if (location != parents_.end()) {
#if DEBUG || DEBUG_RESTRUCTURE
      std::cout << "del_entry_arc: trying to remove arc " << *location
                << " connecting " << (*location)->orig()->description()
                << " to " << (*location)->dest()->description() << endl;
#endif
      parents_.erase(location);
#if DEBUG || DEBUG_RESTRUCTURE
      std::cout << "del_entry_arc: remove arc successful" << endl;
#endif
    } else
      throw std::runtime_error(
          "DGVertex::del_entry_arc() -- the arc doesn't exist");
  } else
    throw std::runtime_error("DGVertex::del_entry_arc() -- no arcs to delete");
}

void DGVertex::detach() {
  // If there are no entry arcs -- then other vertices do not depend on this guy
  // Can safely remove exit arcs
  const unsigned int narcs = num_entry_arcs();
  if (narcs == 0)
    DGVertex::del_exit_arcs();
  else
    throw CannotPerformOperation(
        "DGVertex::detach() -- cannot detach a vertex if it has entry arcs");
}

void DGVertex::prepare_to_traverse() {
  // can_add_arcs_ = false;
  num_tagged_arcs_ = 0;
  scheduled_ = false;
}

unsigned int DGVertex::tag() { return ++num_tagged_arcs_; }

unsigned int DGVertex::num_entry_arcs() const { return parents_.size(); }

unsigned int DGVertex::num_exit_arcs() const { return children_.size(); }

namespace {
struct __ArcDestEqual {
  __ArcDestEqual(const std::shared_ptr<DGVertex>& v) : v_(v) {}
  bool operator()(const std::shared_ptr<DGArc>& a) { return a->dest() == v_; }
  const std::shared_ptr<DGVertex>& v_;
};
}  // namespace

const std::shared_ptr<DGArc>& DGVertex::exit_arc(
    const std::shared_ptr<DGVertex>& v) const {
  static std::shared_ptr<DGArc> nullptr_;
  __ArcDestEqual predicate(v);
  const ArcSetType::const_iterator end = children_.end();
  const ArcSetType::const_iterator pos =
      find_if(children_.begin(), children_.end(), predicate);
  if (pos != end)
    return *pos;
  else
    return nullptr_;
}

void DGVertex::reset() {
  dg_ = 0;
  subtree_ = std::shared_ptr<DRTree>();

  typedef ArcSetType::const_iterator citer;
  typedef ArcSetType::iterator iter;
  const citer end = children_.end();
  for (iter a = children_.begin(); a != end; ++a) {
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
  refs_.resize(0);
}

const std::string& DGVertex::graph_label() const {
  if (!graph_label_.empty())
    return graph_label_;
  else
    throw GraphLabelNotSet("DGVertex::graph_label() -- graph label not set");
}

void DGVertex::set_graph_label(const std::string& label) {
  graph_label_ = label;
}

void DGVertex::refer_this_to(const std::shared_ptr<DGVertex>& V) {
  if (referred_vertex_ != 0) {
    if (referred_vertex_->equiv(V))
      return;
    else
      throw std::logic_error(
          "DGVertex::refer_this_to() -- already referring to some other "
          "vertex");
  }
#if DEBUG
  cout << "DGVertex::refer_this_to() -- vertex " << description()
       << " will refer to " << V->description() << endl;
#endif
  // transfer symbols and addresses to the referred-to index
  if (this->symbol_set() && !V->symbol_set()) V->set_symbol(symbol_);
  if (this->address_set() && !V->address_set()) V->set_symbol(symbol_);
  referred_vertex_ = V.get();
  V->register_reference(this);
}

void DGVertex::register_reference(const DGVertex* referrer) {
  const bool is_new_referrer =
      (std::find(refs_.begin(), refs_.end(), referrer) == refs_.end());
  if (is_new_referrer) {
#if DEBUG
    std::cout << "DGVertex::register_reference() : " << this->description()
              << " has " << refs_.size()
              << " referrers and added new one: " << referrer->description()
              << std::endl;
#endif
    refs_.push_back(referrer);
  } else {
#if DEBUG
    std::cout << "DGVertex::register_reference() : " << this->description()
              << " already has this referrer : " << referrer->description()
              << std::endl;
#endif
  }
  assert(refs_.size() <= 1);
}

const std::string& DGVertex::symbol() const {
  if (referred_vertex_ && referred_vertex_->symbol_set())
    return referred_vertex_->symbol();
  else {
    if (!symbol_.empty())
      return symbol_;
    else {
#if DEBUG
      cout << "DGVertex::symbol() -- symbol not set for " << description()
           << endl;
      if (referred_vertex_)
        cout << "DGVertex::symbol() -- referred_vertex_ = "
             << referred_vertex_->description() << endl;
#endif
      throw SymbolNotSet("DGVertex::symbol() -- symbol not set");
    }
  }
}

bool DGVertex::symbol_set() const {
  if (referred_vertex_)
    return referred_vertex_->symbol_set();
  else
    return !symbol_.empty();
}

void DGVertex::set_symbol(const std::string& symbol) {
  if (referred_vertex_ && referred_vertex_->symbol_set())
    ;  // assert(referred_vertex_->symbol() == symbol);
  else {
    // assert(symbol_.empty());  // should not need to assign twice ... why need
    // to overwrite symbol?
    symbol_ = symbol;
#if DEBUG
    cout << "Set symbol for " << description() << " to " << symbol << endl;
#endif
  }
}

void DGVertex::reset_symbol() { symbol_.clear(); }

DGVertex::Address DGVertex::address() const {
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

bool DGVertex::address_set() const {
  if (referred_vertex_)
    return referred_vertex_->address_set();
  else
    return address_ >= 0;
}

void DGVertex::set_address(const Address& address) { address_ = address; }

void DGVertex::need_to_compute(bool ntc) { need_to_compute_ = ntc; }

bool DGVertex::need_to_compute() const {
  if (referred_vertex_)
    return referred_vertex_->need_to_compute();
  else
    return need_to_compute_;
}

bool DGVertex::precomputed() const {
  if (referred_vertex_)
    return referred_vertex_->precomputed();
  else {
    return this_precomputed();
  }
}

void DGVertex::print(std::ostream& os) const {
  using std::endl;
  std::string prefix("DGVertex::print: ");
  os << prefix << "label = " << label() << endl;
  os << prefix << "this = " << this << endl;
  if (referred_vertex_ != 0) {
    os << prefix << "refers_to = " << referred_vertex_ << endl;
  } else {
    os << prefix << "precomputed = " << precomputed() << endl;
    if (symbol_set()) os << prefix << "symbol = " << symbol() << endl;
    if (address_set()) os << prefix << "address = " << address() << endl;
    os << prefix << "size = " << size() << endl;
    os << prefix << "next to compute = " << postcalc() << endl;
    os << prefix << "nparents = " << num_entry_arcs() << endl;
    unsigned int i = 0;
    for (ArcSetType::const_iterator p = first_entry_arc();
         p != plast_entry_arc(); ++i, ++p)
      os << prefix << "  parent " << i << ": " << (*p)->orig() << endl;
    os << prefix << "nchildren = " << num_exit_arcs() << endl;
    i = 0;
    for (ArcSetType::const_iterator c = first_exit_arc(); c != plast_exit_arc();
         ++i, ++c)
      os << prefix << "  child " << i << ": " << (*c)->dest() << endl;
    os << prefix << "ntags = " << num_tagged_arcs_ << endl;
  }
}

void DGVertex::unregister() const {}

////

bool UnrolledIntegralSet::operator()(const std::shared_ptr<DGVertex>& V) {
  const unsigned int outdegree = V->num_exit_arcs();
  if (outdegree == 0) return false;

  const std::shared_ptr<DGArc> arc0 = *(V->first_exit_arc());
  // Is this DGArcRR?
  const std::shared_ptr<DGArcRR> arcrr =
      std::dynamic_pointer_cast<DGArcRR, DGArc>(arc0);
  if (arcrr == 0) return false;
  // Is this DGArcRR<IntegralSet_to_Integral>? If invariant_type() is false,
  // then yes
  return !arcrr->rr()->invariant_type();
}

bool NotUnrolledIntegralSet::operator()(const std::shared_ptr<DGVertex>& V) {
  return !UnrolledIntegralSet()(V);
}

bool IntegralInTargetIntegralSet::operator()(
    const std::shared_ptr<DGVertex>& V) {
  const unsigned int indegree = V->num_entry_arcs();
  if (indegree != 1) return false;
  auto parent = (*(V->first_entry_arc()))->orig();
  if (parent->is_a_target()) return UnrolledIntegralSet()(parent);
  return false;
}
