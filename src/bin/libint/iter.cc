
#include <iter.h>

using namespace std;
using namespace libint2;

//
// SubIterator
//
SubIterator::SubIterator()
{
}

SubIterator::~SubIterator()
{
}

#ifdef ALLOW_PARTIALLY_SPECIALIZED_NESTED_TEMPLATES

//
// SubIteratorBase::init_subobj<CGF>
//
template <class T, class P>
template <>
void
SubIteratorBase<T,P>::init_subobj<CGF>()
{
  subobj_.push_back(obj_);
}

template <class T, class P>
template <>
void
SubIteratorBase<T,P>::delete_subobj<CGF>()
{
}

//
// SubIteratorBase::init_subobj<CGShell>
//
template <class T, class P>
template <>
void
SubIteratorBase<T,P>::init_subobj<CGShell>()
{
  P::cgshell_to_cgfvector(obj_,subobj_);
}

template <class T, class P>
template <>
void
SubIteratorBase<T,P>::delete_subobj<CGShell>()
{
  int nelem = subobj_.size();
  for(int i=0; i<nelem; i++)
    subobj_[i]->~CGF();
}

#endif

