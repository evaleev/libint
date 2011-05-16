
#include <cstdlib>
#include <ctime>
#include <dg.h>
#include <dg.templ.h>
#include <tactic.h>

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
      else {
	std::cout << "ZeroNewVerticesTactic::optimal_rr: not optimal: " << stack[i]->label() << std::endl;
	SafePtr<DGVertex> target = stack[i]->rr_target();
	const unsigned int nchildren = stack[i]->num_children();
	for(unsigned int c=0; c<nchildren; ++c) {
	  SafePtr<DGVertex> child = stack[i]->rr_child(c);
	  std::cout << "  child " << c << ": " << child->label() << std::endl;
	}
      }
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
    const unsigned long range = 1ul<<32 - 1;
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


