
#include <rr.h>
#include <traits.h>

using namespace std;
using namespace libint2;

//
// StdLibintTraits<CGShell>
//
template <>
void
StdLibintTraits<CGShell>::init_subobj(const CGShell* cgshell, vector<const CGF*>& cgfs)
{
  unsigned int am = cgshell->qn();
  unsigned int qn[3];
  for(unsigned int i=0; i<=am; i++) {
    qn[0] = am - i;
    for(unsigned int j=0; j<=i; j++) {
      qn[1] = i - j;
      qn[2] = j;

      cgfs.push_back(new CGF(qn));
    }
  }
}

template <>
void
StdLibintTraits<CGShell>::dealloc_subobj(vector<const CGF*>& subobj)
{
  int nelem = subobj.size();
  for(int i=0; i<nelem; i++)
    subobj[i]->~CGF();
}

//
// StdLibintTraits<CGF>
//
template <>
void
StdLibintTraits<CGF>::init_subobj(const CGF* obj, vector<const CGF*>& subobj)
{
  subobj.push_back(obj);
}

template <>
void
StdLibintTraits<CGF>::dealloc_subobj(vector<const CGF*>& subobj)
{
}

//
// StdLibintTraits<TwoERep>
//
/*
template <>
void
StdLibintTraits<TwoERep>::init_subobj(const TwoERep* obj, vector<const TwoERep*>& subobj)
{
  subobj.push_back(obj);
}

template <>
void
StdLibintTraits<TwoERep>::dealloc_subobj(vector<const TwoERep*>& subobj)
{
}
*/
