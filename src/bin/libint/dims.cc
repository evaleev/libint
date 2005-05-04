
#include <default_params.h>
#include <dims.h>

using namespace std;
using namespace libint2;

SafePtr<ImplicitDimensions>
ImplicitDimensions::default_dims_(new ImplicitDimensions(1,1,StaticDefinitions::default_vector_length));

SafePtr<ImplicitDimensions>
ImplicitDimensions::default_dims()
{
  return default_dims_;
}


ImplicitDimensions::ImplicitDimensions(const SafePtr<Entity>& high,
                                       const SafePtr<Entity>& low,
                                       const SafePtr<Entity>& vecdim) :
  high_(high), low_(low), vecdim_(vecdim)
  {
    init_();
  }

ImplicitDimensions::ImplicitDimensions() :
  high_(SafePtr<Entity>(new RTimeEntity<EntityTypes::Int>("highdim"))),
  low_(SafePtr<Entity>(new RTimeEntity<EntityTypes::Int>("lowdim"))),
  vecdim_(SafePtr<Entity>(new RTimeEntity<EntityTypes::Int>("libint->veclen")))
  {
    init_();
  }

ImplicitDimensions::ImplicitDimensions(int high, int low, int vec) :
  high_(SafePtr<Entity>(new CTimeEntity<int>("highdim",high))),
  low_(SafePtr<Entity>(new CTimeEntity<int>("lowdim",low))),
  vecdim_(SafePtr<Entity>(new CTimeEntity<int>("vecdim",vec)))
  {
    init_();
  }

void
ImplicitDimensions::init_()
{
  SafePtr< CTimeEntity<int> > cptr = dynamic_pointer_cast<CTimeEntity<int>,Entity>(high_);
  if (cptr != 0)
    high_is_static_ = true;
  else
    high_is_static_ = false;
  cptr = dynamic_pointer_cast<CTimeEntity<int>,Entity>(low_);
  if (cptr != 0)
    low_is_static_ = true;
  else
    low_is_static_ = false;
  cptr = dynamic_pointer_cast<CTimeEntity<int>,Entity>(vecdim_);
  if (cptr != 0)
    vecdim_is_static_ = true;
  else
    vecdim_is_static_ = false;
}

