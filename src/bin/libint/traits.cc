
#include <traits.h>
#include <rr.h>

using namespace std;
using namespace libint2;

template <>
void
StdLibintTraits<CGShell>::init_subobj(const CGShell* cgshell, vector<const CGF*>& cgfs)
{
  unsigned int am = cgshell->qn();
  unsigned int qn[3];
  for(int i=0; i<=am; i++) {
    qn[0] = am - i;
    for(int j=0; j<=i; j++) {
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

StdLibintTraits< VectorBraket<CGShell> > a;

