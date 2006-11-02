
#include <graph_registry.h>

using namespace libint2;

GraphRegistry::GraphRegistry() :
  accumulate_targets_(false), return_targets_(true), can_unroll_(true), do_cse_(false),
  stack_name_("inteval->stack")
{
}

GraphRegistry::~GraphRegistry()
{
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
