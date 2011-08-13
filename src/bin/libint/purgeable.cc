
#include <purgeable.h>
#include <algorithm>
#include <cassert>

using namespace libint2;

PurgeableStacks*
PurgeableStacks::instance_ = 0;

PurgeableStacks*
PurgeableStacks::Instance() {
  if (instance_ == 0)
    instance_ = new PurgeableStacks;
  return instance_;
}

void
PurgeableStacks::register_stack(stack_type* stack) {
  typedef std::vector<stack_type*>::const_iterator citer;
  citer v = std::find(stacks_.begin(),stacks_.end(),stack);
  assert(v == stacks_.end());
  stacks_.push_back(stack);
}

void
PurgeableStacks::purge() {
  typedef std::vector<stack_type*>::iterator iter;
  for(iter s=stacks_.begin(); s!=stacks_.end(); ++s)
    (*s)->purge();
}
