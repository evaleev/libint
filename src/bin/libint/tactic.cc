/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

#include <cstdlib>
#include <ctime>
#include <dg.h>
#include <tactic.h>
#include <rr.h>
#include <master_ints_list.h>

using namespace std;
using namespace libint2;

/*
FirstChoiceTactic::RR
FirstChoiceTactic::optimal_rr(const rr_stack& stack) const {
  if (!stack.empty())
    return stack[0];
  else
    return RR();
}
*/

FewestNewVerticesTactic::RR
FewestNewVerticesTactic::optimal_rr(const rr_stack& stack) const {
  if (!stack.empty()) {
    unsigned int best_rr = 0;
    unsigned int min_nchildren = 1000000000;
    // Loop over all RRs and find the one with the fewest children
    for(unsigned int i=0; i<stack.size(); i++) {
      unsigned int nchildren = dg_->num_children_on(stack[i]);
      if (nchildren < min_nchildren) {
        min_nchildren = nchildren;
        best_rr = i;
      }
    }
    return stack[best_rr];
  }
  else
    // Else return null pointer
  return RR();
}

ZeroNewVerticesTactic::RR
ZeroNewVerticesTactic::optimal_rr(const rr_stack& stack) const {
  if (!stack.empty()) {
    // Loop over all RRs and find the one with zero children
    for(unsigned int i=0; i<stack.size(); i++) {
      const RR& rr = stack[i];
      const unsigned int nchildren = rr->num_children();
      unsigned int nchildren_on_dg = dg_->num_children_on(rr);
      if (nchildren == nchildren_on_dg) {
        return rr;
      }
//      else {
//        std::cout << "ZeroNewVerticesTactic::optimal_rr: not optimal: " << stack[i]->label() << std::endl;
//        SafePtr<DGVertex> target = stack[i]->rr_target();
//        const unsigned int nchildren = stack[i]->num_children();
//        for(unsigned int c=0; c<nchildren; ++c) {
//          SafePtr<DGVertex> child = stack[i]->rr_child(c);
//          std::cout << "  child " << c << ": " << child->label() << std::endl;
//        }
//
//        {
//          VertexPrinter vp(std::cout);
//          dg_->foreach(vp);
//        }
//
//      }
    }
    throw std::logic_error("ZeroNewVerticesTactic -- no RRs found that add zero new vertices. Probably used by mistake");
  }
  else
    // Else return null pointer
  return RR();
}

RandomChoiceTactic::RandomChoiceTactic() : Tactic()
{
  // Initialize state randomly
  time_t crap;
  srandom(time(&crap));
}

RandomChoiceTactic::RR
RandomChoiceTactic::optimal_rr(const rr_stack& stack) const {
  if (!stack.empty()) {
    unsigned int size = stack.size();
    unsigned long rand = random();
    const unsigned long range = RAND_MAX;
    long choice = (long)(rand * size - 1)/range;
    return stack[choice];
  }
  else
    return RR();
}

NullTactic::RR
NullTactic::optimal_rr(const rr_stack& stack) const {
  return RR();
}

ParticleDirectionTactic::RR
ParticleDirectionTactic::optimal_rr(const rr_stack& stack) const {
  //  std::cout << "in ParticleDirectionTactic::optimal_rr : increase_=" << increase_ << std::endl;

  // try to find the first RR with matching direction
  for (auto& t : stack) {
//    std::cout << " rr=" << t->label() << std::endl;
    if (t->partindex_direction() == +1 && increase_)
      return t;
    if (t->partindex_direction() == -1 && not increase_)
      return t;
//    std::cout << "**not selected**" << std::endl;
  }

  // if failed to find an RR with matching direction, choose the first non-directional
  for (auto& t : stack) {
    if (t->partindex_direction() == 0)
      return t;
  }

  // if all failed, return an empty RR
  return RR();
}

FourCenter_OS_Tactic::RR
FourCenter_OS_Tactic::optimal_rr(const rr_stack& stack) const {

  if (stack.empty())
    return RR();

  // grab the quantum numbers of the target set of these RRs
  unsigned lbra0, lket0, lbra1, lket1;
  SafePtr<TwoPRep_11_11_sq> abcd_ptr = dynamic_pointer_cast<TwoPRep_11_11_sq>(stack[0]->rr_target());
  if (abcd_ptr) {
    lbra0 = abcd_ptr->bra(0,0).norm();
    lbra1 = abcd_ptr->bra(1,0).norm();
    lket0 = abcd_ptr->ket(0,0).norm();
    lket1 = abcd_ptr->ket(1,0).norm();
  }
  else {
    SafePtr<TwoPRep_11_11_int> abcd_ptr = dynamic_pointer_cast<TwoPRep_11_11_int>(stack[0]->rr_target());
    if (abcd_ptr) {
      lbra0 = abcd_ptr->bra(0,0).norm();
      lbra1 = abcd_ptr->bra(1,0).norm();
      lket0 = abcd_ptr->ket(0,0).norm();
      lket1 = abcd_ptr->ket(1,0).norm();
    }
    else {
      assert(false); // should not be possible
    }
  }

  auto l0 = lbra0_ + lket0_;
  auto l1 = lbra1_ + lket1_;
  auto l1_ge_l0 = l1 >= l0;
  auto use_itr = l0 > 0 && l1 > 0;

  // try to apply ITR first
  if (use_itr) {
    for (auto& t : stack) {

      if (t->partindex_direction() == 0) // skip all non-ITR RRs
        continue;

      //
      // this ITR is useful if it shifts the quanta towards the particle with greater total quanta
      if (t->partindex_direction() == +1 && l1_ge_l0)
        return t;
      if (t->partindex_direction() == -1 && not l1_ge_l0)
        return t;
    }
  }
  // else use the non-ITR relation with smallest size of children (to reduce the memory bandwidth demand)
  RR result;
  size_t max_result_size = std::numeric_limits<size_t>::max();
  size_t nbests = 0;
  for (auto& t : stack) {
    if (t->partindex_direction() == 0) { // skip all ITR RRs
      if (t->size_of_children() < max_result_size)
        ++nbests;
      if (t->size_of_children() < max_result_size) {
        result = t;
        max_result_size = t->size_of_children();
        nbests = 1;
      }
    }
  }

  if (nbests > 1) {
    std::cout << "FourCenter_OS_Tactic: found more than one RR with same (optimal) size of children" << std::endl;
    assert(nbests == 1);
  }

  return result;
}
