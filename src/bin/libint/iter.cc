
#include <iter.h>

using namespace std;
using namespace libint2;


/** SetIterator<CGShell> is a specialization of SetIterator<>.
 */
template <>
const unsigned int
SetIterator<CGShell>::num_iter() const
{
  return obj_->num_bf();
}

/** SetIterator<CGF> is a specialization of SetIterator<>.
 */
template <>
const unsigned int
SetIterator<CGF>::num_iter() const
{
  return obj_->num_bf();
}

/** SetIterator<CGF> is a specialization of SetIterator<>.
 */
template <>
const SetIterator<CGF>::iter_type*
SetIterator<CGF>::first()
{
  return obj_;
}

/** SetIterator<CGF> is a specialization of SetIterator<>.
 */
template <>
const SetIterator<CGF>::iter_type*
SetIterator<CGF>::next()
{
  return 0;
}
  
