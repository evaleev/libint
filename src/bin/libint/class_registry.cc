
#include <class_registry.h>

using namespace libint2;

ClassRegistry::ClassRegistry() :
  nclasses_(0)
{
}

ClassRegistry*
ClassRegistry::registry_ = 0;

ClassRegistry&
ClassRegistry::Instance()
{
  if (!registry_)
    registry_ = new ClassRegistry;
  return *registry_;
}
