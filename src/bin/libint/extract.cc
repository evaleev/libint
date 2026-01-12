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
#include <entity.h>
#include <extract.h>
#include <intset_to_ints.h>
#include <rr.h>
#include <uncontract.h>

#include <iostream>
#include <map>

using namespace std;
using namespace libint2;

void ExtractExternSymbols::operator()(const VertexPtr& v) {
#if DEBUG
  std::cout << "ExtractExternSymbols::operator() -- v = " << v->description()
            << std::endl;
#endif
  if (v->precomputed()) {
    // discard compile-time entities
    {
      typedef CTimeEntity<double> cdouble;
      std::shared_ptr<cdouble> ptr_cast =
          std::dynamic_pointer_cast<cdouble, DGVertex>(v);
      if (ptr_cast) {
        return;
      }
    }

    // discard unrolled integral sets composed of precomputed integrals
    {
      std::shared_ptr<DGArcRR> arcrr;
      if (v->size() == 1 && v->num_exit_arcs() == 1 &&
          ((arcrr = std::dynamic_pointer_cast<DGArcRR, DGArc>(
                *(v->first_exit_arc()))) != 0
               ? std::dynamic_pointer_cast<IntegralSet_to_Integrals_base,
                                           RecurrenceRelation>(arcrr->rr()) != 0
               : false) &&
          (*(v->first_exit_arc()))->dest()->precomputed()) {
        return;
      }
    }

    map_[v->label()] = true;
  }
}

const ExtractExternSymbols::Symbols& ExtractExternSymbols::symbols() const {
  symbols_.clear();
  typedef LabelMap::const_iterator citer;
  citer end = map_.end();
  for (citer l = map_.begin(); l != end; ++l) {
    symbols_.push_back(l->first);
  }
  // symbols_.sort();
  return symbols_;
}

////

void ExtractRR::operator()(const VertexPtr& v) {
  if (v->num_exit_arcs() != 0) {
    std::shared_ptr<DGArc> arc = *(v->first_exit_arc());
    std::shared_ptr<DGArcRR> arc_rr =
        std::dynamic_pointer_cast<DGArcRR, DGArc>(arc);
    if (arc_rr != 0) {
      std::shared_ptr<RecurrenceRelation> rr = arc_rr->rr();
      std::shared_ptr<IntegralSet_to_Integrals_base> iset_to_i =
          std::dynamic_pointer_cast<IntegralSet_to_Integrals_base,
                                    RecurrenceRelation>(rr);
      std::shared_ptr<Uncontract_Integral_base> unc_i =
          std::dynamic_pointer_cast<Uncontract_Integral_base,
                                    RecurrenceRelation>(rr);
      if (iset_to_i == 0 && unc_i == 0) {
        const std::shared_ptr<RRStack>& rrstack = RRStack::Instance();
        // RRStack must be guaranteed to have this rr
        const RRStack::value_type rrstackvalue = rrstack->find(rr);
        const RRid rrid = rrstackvalue.first;
        map_[rrid] = true;
      }
    }
  }
}

const ExtractRR::RRList& ExtractRR::rrlist() const {
  rrlist_.clear();
  typedef RRMap::const_iterator citer;
  citer end = map_.end();
  for (citer rr = map_.begin(); rr != end; ++rr) {
    rrlist_.push_back(rr->first);
  }
  // rrlist_.sort();
  return rrlist_;
}
