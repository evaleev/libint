
#include <map>
#include <iostream>
#include <dgvertex.h>
#include <entity.h>
#include <rr.h>
#include <intset_to_ints.h>
#include <extract.h>

using namespace std;
using namespace libint2;

void
ExtractExternSymbols::operator()(const VertexPtr& v)
{
#if DEBUG
  std::cout << "ExtractExternSymbols::operator() -- v = " << v->description() << std::endl;
#endif
  if (v->precomputed()) {

    // discard compile-time entities
    {
      typedef CTimeEntity<double> cdouble;
      SafePtr<cdouble> ptr_cast = dynamic_pointer_cast<cdouble,DGVertex>(v);
      if (ptr_cast) {
        return;
      }
    }
    
    // discard unrolled integral sets composed of precomputed integrals
    {
      SafePtr<DGArcRR> arcrr;
      if (v->size() == 1 && v->num_exit_arcs() == 1 &&
          ( (arcrr = dynamic_pointer_cast<DGArcRR,DGArc>(*(v->first_exit_arc()))) != 0 ?
              dynamic_pointer_cast<IntegralSet_to_Integrals_base,RecurrenceRelation>(arcrr->rr()) != 0 :
              false ) &&
          (*(v->first_exit_arc()))->dest()->precomputed()
         ) {
        return;
      }
    }

    map_[v->label()] = true;
  }

}

const ExtractExternSymbols::Symbols&
ExtractExternSymbols::symbols() const
{
  symbols_.clear();
  typedef LabelMap::const_iterator citer;
  citer end = map_.end();
  for(citer l=map_.begin(); l!=end; ++l) {
    symbols_.push_back(l->first);
  }
  //symbols_.sort();
  return symbols_;
}

////

void
ExtractRR::operator()(const VertexPtr& v)
{
  if (v->num_exit_arcs() != 0) {
    SafePtr<DGArc> arc = *(v->first_exit_arc());
    SafePtr<DGArcRR> arc_rr = dynamic_pointer_cast<DGArcRR,DGArc>(arc);
    if (arc_rr != 0) {
      SafePtr<RecurrenceRelation> rr = arc_rr->rr();
      SafePtr<IntegralSet_to_Integrals_base> iset_to_i = dynamic_pointer_cast<IntegralSet_to_Integrals_base,RecurrenceRelation>(rr);
      if (iset_to_i == 0) {
	const SafePtr<RRStack>& rrstack = RRStack::Instance();
	// RRStack must be guaranteed to have this rr
	const RRStack::value_type rrstackvalue = rrstack->find(rr);
	const RRid rrid = rrstackvalue.first;
	map_[rrid] = true;
      }
    }
  }

}

const ExtractRR::RRList&
ExtractRR::rrlist() const
{
  rrlist_.clear();
  typedef RRMap::const_iterator citer;
  citer end = map_.end();
  for(citer rr=map_.begin(); rr!=end; ++rr) {
    rrlist_.push_back(rr->first);
  }
  //rrlist_.sort();
  return rrlist_;
}
