
#include <graph_registry.h>

using namespace libint2;

GraphRegistry::GraphRegistry() :
  accumulate_targets_(false), return_targets_(true), unroll_threshold_(1), uncontract_(false), ignore_missing_prereqs_(false),
  do_cse_(false), condense_expr_(false), stack_name_("inteval->stack")
{
}

GraphRegistry::~GraphRegistry()
{
}

GraphRegistry*
GraphRegistry::clone() const {
  GraphRegistry* gr = new GraphRegistry;
  *gr = *this;
  return gr;
}

////

InternalGraphRegistry::InternalGraphRegistry() :
  accumulate_targets_directly_(false),
  size_of_target_accum_(0)
{
}

InternalGraphRegistry::~InternalGraphRegistry()
{
}
