
#include <policy.h>
#include <rr.h>

using namespace std;
using namespace libint2;

/*
 Definition of a generic StdLibintTDPolicy is provided in policy.h
 */

/**
StdLibintTDPolicy<CGShell>::init_subobj initializes CGFs in canonical order.
 The functions in order are produced using the following C++ loop:

 for(int i=0; i<=am; i++) {
   qn[0] = am - i;
   for(int j=0; j<=i; j++) {
     qn[1] = i - j;
     qn[2] = j;
   }
 }

 where am is the angular momentum of the shell and qn[3] are the x, y, and z
 exponents.
 */

template <>
void
StdLibintTDPolicy<CGShell>::init_subobj(const StdLibintTDPolicy<CGShell>::obj_stype& cgshell,
vector<StdLibintTDPolicy<CGShell>::subobj_stype>& cgfs)
{
  unsigned int am = TypeTraits<CGShell>::const_ref(cgshell).qn();
  unsigned int qn[3];
  for(unsigned int i=0; i<=am; i++) {
    qn[0] = am - i;
    for(unsigned int j=0; j<=i; j++) {
      qn[1] = i - j;
      qn[2] = j;
      
      subobj_stype cgf(qn);
      cgfs.push_back(cgf);
    }
  }
}

template <>
void
StdLibintTDPolicy<CGShell>::dealloc_subobj(vector<StdLibintTDPolicy<CGShell>::subobj_stype>& subobj)
{
}

