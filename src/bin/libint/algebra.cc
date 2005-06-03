
#include <algebra.h>

using namespace libint2;

template <>
void
AlgebraicOperator<DGVertex>::add_exit_arc(const SafePtr<DGArc>& a)
{
  DGVertex::add_exit_arc(a);
  if (left_->equiv(a->dest()))
    left_ = a->dest();
  else if (right_->equiv(a->dest()))
    right_ = a->dest();
  else
    throw std::runtime_error("AlgebraicOperator<DGVertex>::add_exit_arc -- trying to add an arc to a vertex not equivalent to either argument.");
}
