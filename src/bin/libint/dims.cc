
#include <dims.h>

using namespace std;
using namespace libint2;

SafePtr<ImplicitDimensions>
ImplicitDimensions::default_dims_(new ImplicitDimensions(1,1));

SafePtr<ImplicitDimensions>
ImplicitDimensions::default_dims()
{
  return default_dims_;
}

void
ImplicitDimensions::init_()
{
  SafePtr< CTimeEntity<int> > cptr = boost::shared_ptr< CTimeEntity<int> >(high_,boost::detail::dynamic_cast_tag());
  if (cptr != 0)
    high_is_static_ = true;
  else
    high_is_static_ = false;
  cptr = boost::shared_ptr< CTimeEntity<int> >(low_,boost::detail::dynamic_cast_tag());
  if (cptr != 0)
    low_is_static_ = true;
  else
    low_is_static_ = false;
}

