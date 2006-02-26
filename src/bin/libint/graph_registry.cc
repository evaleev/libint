
#include <graph_registry.h>

using namespace libint2;

GraphRegistry::GraphRegistry() :
  can_unroll_(true), do_cse_(false),
  stack_name_("libint->stack")
{
}

GraphRegistry::~GraphRegistry()
{
}

