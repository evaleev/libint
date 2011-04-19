
#include <policy.h>
#include <rr.h>
#include <cgshell_ordering.h>
#include <cgshellinfo.h>

using namespace std;

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

namespace {
  /// returns an xyz which is not a and not b
  int notxyz(int a, int b);
  /// returns an ordered pair of xyz (i.e. xy, not yx) which is not a
  std::pair<int,int> notxyz(int a);
}

namespace libint2 {

template <>
void
StdLibintTDPolicy<CGShell>::init_subobj(const StdLibintTDPolicy<CGShell>::obj_stype& cgshell,
vector<StdLibintTDPolicy<CGShell>::subobj_stype>& cgfs)
{
  unsigned int am = TypeTraits<CGShell>::const_ref(cgshell).qn();
  unsigned int qn[3] = {0, 0, 0};
  int lx, ly, lz;
  FOR_CART(lx,ly,lz,am)
    qn[0] = lx;
    qn[1] = ly;
    qn[2] = lz;
    subobj_stype cgf(qn);
    cgf.deriv() = cgshell.deriv();
    cgfs.push_back(cgf);
  END_FOR_CART
}

template <>
void
StdLibintTDPolicy<CGShell>::dealloc_subobj(vector<StdLibintTDPolicy<CGShell>::subobj_stype>& subobj)
{
}

};

////

namespace {
  int notxyz(int a, int b) {
    if (a == b)
      throw libint2::ProgrammingError("notxyz(a,b) -- a equals b");
    int amax = std::max(a,b);
    int amin = std::min(a,b);
    if (amin == 0 && amax == 1)
      return 2;
    if (amin == 0 && amax == 2)
      return 1;
    if (amin == 1 && amax == 2)
      return 0;
    abort(); // unreachable
  }

  std::pair<int,int> notxyz(int a) {
    switch(a) {
    case 0: return make_pair(1,2); break;
    case 1: return make_pair(0,2); break;
    case 2: return make_pair(0,1); break;
    }
    abort(); // unreachable
  }
}
