
#include <algebra.h>
#include <global_macros.h>

namespace libint2 {

  template <>
  void
  AlgebraicOperator<DGVertex>::add_exit_arc(const SafePtr<DGArc>& a)
  {
    DGVertex::add_exit_arc(a);
#if CHECK_SAFETY
    if (num_exit_arcs() > 2)
      throw std::runtime_error("AlgebraicOperator<DGVertex>::add_exit_arc() -- number of exit arcs is now greater than 2!");
#endif
    if (left_->equiv(a->dest()))
      left_ = a->dest();
    else if (right_->equiv(a->dest()))
      right_ = a->dest();
    else
      throw std::runtime_error("AlgebraicOperator<DGVertex>::add_exit_arc -- trying to add an arc to a vertex not equivalent to either argument.");
  }
}

