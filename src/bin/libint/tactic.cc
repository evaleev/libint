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

#include <dg.h>
#include <master_ints_list.h>
#include <rr.h>
#include <tactic.h>

#include <cassert>
#include <cstdlib>
#include <ctime>

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

FewestNewVerticesTactic::RR FewestNewVerticesTactic::optimal_rr(
    const rr_stack& stack) const {
  if (!stack.empty()) {
    unsigned int best_rr = 0;
    unsigned int min_nchildren = 1000000000;
    // Loop over all RRs and find the one with the fewest new children
    for (unsigned int i = 0; i < stack.size(); i++) {
      int nchildren = stack[i]->num_children();
      int nchildren_on_dg = dg_->num_children_on(stack[i]);
      assert(nchildren >= nchildren);
      int nchildren_new = nchildren - nchildren_on_dg;
      if (nchildren_new < min_nchildren) {
        min_nchildren = nchildren_new;
        best_rr = i;
      }
    }
    return stack[best_rr];
  } else
    // Else return null pointer
    return RR();
}

ZeroNewVerticesTactic::RR ZeroNewVerticesTactic::optimal_rr(
    const rr_stack& stack) const {
  if (!stack.empty()) {
    // Loop over all RRs and find the first one with zero children
    for (unsigned int i = 0; i < stack.size(); i++) {
      const RR& rr = stack[i];
      const unsigned int nchildren = rr->num_children();
      unsigned int nchildren_on_dg = dg_->num_children_on(rr);
      if (nchildren == nchildren_on_dg) {
        return rr;
      }
      //      else {
      //        std::cout << "ZeroNewVerticesTactic::optimal_rr: not optimal: "
      //        << stack[i]->label() << std::endl; std::shared_ptr<DGVertex>
      //        target = stack[i]->rr_target(); const unsigned int nchildren =
      //        stack[i]->num_children(); for(unsigned int c=0; c<nchildren;
      //        ++c) {
      //          std::shared_ptr<DGVertex> child = stack[i]->rr_child(c);
      //          std::cout << "  child " << c << ": " << child->label() <<
      //          std::endl;
      //        }
      //
      //        {
      //          VertexPrinter vp(std::cout);
      //          dg_->foreach(vp);
      //        }
      //
      //      }
    }
    throw std::logic_error(
        "ZeroNewVerticesTactic -- no RRs found that add zero new vertices. "
        "Probably used by mistake");
  } else
    // Else return null pointer
    return RR();
}

RandomChoiceTactic::RandomChoiceTactic() : Tactic() {
  // Initialize state randomly
  time_t crap;
  srandom(time(&crap));
}

RandomChoiceTactic::RR RandomChoiceTactic::optimal_rr(
    const rr_stack& stack) const {
  if (!stack.empty()) {
    unsigned int size = stack.size();
    unsigned long rand = random();
    const unsigned long range = RAND_MAX;
    long choice = (long)(rand * size - 1) / range;
    return stack[choice];
  } else
    return RR();
}

NullTactic::RR NullTactic::optimal_rr(const rr_stack& stack) const {
  return RR();
}

ParticleDirectionTactic::RR ParticleDirectionTactic::optimal_rr(
    const rr_stack& stack) const {
  //  std::cout << "in ParticleDirectionTactic::optimal_rr : increase_=" <<
  //  increase_ << std::endl;

  // try to find the first RR with matching direction
  for (auto& t : stack) {
    //    std::cout << " rr=" << t->label() << std::endl;
    if (t->partindex_direction() == +1 && increase_) return t;
    if (t->partindex_direction() == -1 && not increase_) return t;
    //    std::cout << "**not selected**" << std::endl;
  }

  // if failed to find an RR with matching direction, choose the first
  // non-directional
  for (auto& t : stack) {
    if (t->partindex_direction() == 0) return t;
  }

  // if all failed, return an empty RR
  return RR();
}

TwoCenter_OS_Tactic::RR TwoCenter_OS_Tactic::optimal_rr(
    const rr_stack& stack) const {
  RR result;

  if (!stack.empty()) {
    auto use_hrr = lbra0_ > 0 && lket0_ > 0;

    // try to apply HRR first
    if (use_hrr) {
      for (auto& t : stack) {
        if (t->braket_direction() ==
            BraketDirection::None)  // skip all non-HRR RRs
          continue;

        if (class_debug()) {
          std::cout << "TwoCenter_OS_Tactic: considering " << t->label()
                    << std::endl;
        }

        //
        // this HRR is useful if it shifts the quanta towards the particle with
        // greater total quanta
        if (t->braket_direction() == BraketDirection::KetToBra &&
            lbra0_ >= lket0_) {
          result = t;
          continue;
        } else if (t->braket_direction() == BraketDirection::BraToKet &&
                   lbra0_ < lket0_) {
          result = t;
          continue;
        }
      }
    }

    if (!result) {  // else use the non-HRR relation with smallest size of
                    // children
                    // (to reduce the memory bandwidth demand)
      size_t max_result_size = std::numeric_limits<size_t>::max();
      LIBINT_MAYBE_UNUSED size_t nties = 0;
      for (auto& t : stack) {
        if (t->braket_direction() ==
            BraketDirection::None) {  // skip all HRR RRs
          if (class_debug()) {
            std::cout << "TwoCenter_OS_Tactic: considering " << t->label()
                      << std::endl;
          }
          if (t->size_of_children() == max_result_size)
            ++nties;
          else if (t->size_of_children() < max_result_size) {
            result = t;
            max_result_size = t->size_of_children();
            nties = 0;
          }
        }
      }

      // TODO determine how to resolve ties in 2-center OS tactic
      //      if (nties > 1) {
      //        std::cout << "TwoCenter_OS_Tactic: found more than one RR with
      //        same "
      //                     "(optimal) size of children"
      //                  << std::endl;
      //        assert(nties == 0);
      //      }
    }

    if (class_debug()) {
      if (result)
        std::cout << "TwoCenter_OS_Tactic: picked " << result->label()
                  << std::endl;
      else
        std::cout << "TwoCenter_OS_Tactic: picked none" << std::endl;
    }
  }

  return result;
}

FourCenter_OS_Tactic::RR FourCenter_OS_Tactic::optimal_rr(
    const rr_stack& stack) const {
  RR result;

  if (!stack.empty()) {
    auto l0 = lbra0_ + lket0_;
    auto l1 = lbra1_ + lket1_;
    auto l1_ge_l0 = l1 >= l0;
    auto use_itr = l0 > 0 && l1 > 0;

    // try to apply ITR first
    if (use_itr) {
      for (auto& t : stack) {
        if (t->partindex_direction() == 0)  // skip all non-ITR RRs
          continue;

        if (class_debug()) {
          std::cout << "FourCenter_OS_Tactic: considering " << t->label()
                    << std::endl;
        }

        //
        // this ITR is useful if it shifts the quanta towards the particle with
        // greater total quanta
        if (t->partindex_direction() == +1 && l1_ge_l0) {
          result = t;
          continue;
        } else if (t->partindex_direction() == -1 && not l1_ge_l0) {
          result = t;
          continue;
        }
      }
    }

    if (!result) {
      // else use the non-ITR relation with smallest size of children (to reduce
      // the memory bandwidth demand) note that there is no need to check if
      // transfer direction matches the strategic direction since only non-ITR
      // 2-body OS strategies will include transfers in single direction
      size_t max_result_size = std::numeric_limits<size_t>::max();
      LIBINT_MAYBE_UNUSED size_t nties = 0;
      for (auto& t : stack) {
        if (t->partindex_direction() == 0) {  // skip all ITR RRs

          if (class_debug()) {
            std::cout << "FourCenter_OS_Tactic: considering " << t->label()
                      << std::endl;
          }

          if (t->size_of_children() == max_result_size)
            ++nties;
          else if (t->size_of_children() < max_result_size) {
            result = t;
            max_result_size = t->size_of_children();
            nties = 0;
          }
        }
      }

      // TODO determine how to resolve ties in 4-center OS tactic
      //      if (nties > 0) {
      //        std::cout << "FourCenter_OS_Tactic: found more than one RR with
      //        same "
      //                     "(optimal) size of children"
      //                  << std::endl;
      //        assert(nties == 0);
      //      }
    }

    if (class_debug()) {
      if (result)
        std::cout << "FourCenter_OS_Tactic: picked " << result->label()
                  << std::endl;
      else
        std::cout << "FourCenter_OS_Tactic: picked none" << std::endl;
    }
  }

  return result;
}
