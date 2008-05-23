
#include <policy.h>
#include <rr.h>
#include <cgshell_ordering.h>

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

#if 0
#if LIBINT_CGSHELL_ORDERING == LIBINT__CGSHELL__ORDERING_STANDARD
template <>
void
StdLibintTDPolicy<CGShell>::init_subobj(const StdLibintTDPolicy<CGShell>::obj_stype& cgshell,
vector<StdLibintTDPolicy<CGShell>::subobj_stype>& cgfs)
{
  unsigned int am = TypeTraits<CGShell>::const_ref(cgshell).qn();
  unsigned int qn[3] = {0, 0, 0};
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
#endif
#if LIBINT_CGSHELL_ORDERING == LIBINT__CGSHELL__ORDERING_INTV3
template <>
void
StdLibintTDPolicy<CGShell>::init_subobj(const StdLibintTDPolicy<CGShell>::obj_stype& cgshell,
vector<StdLibintTDPolicy<CGShell>::subobj_stype>& cgfs)
{
  const unsigned int am = TypeTraits<CGShell>::const_ref(cgshell).qn();
  unsigned int qn[3] = {0, 0, 0};
  qn[1] = am;
  const unsigned int nbf = cgshell.num_bf();
  for(unsigned int bf=0; bf < nbf; ++bf) {
    subobj_stype cgf(qn);
    cgfs.push_back(cgf);
    
    if (qn[2] < am - qn[0])
      ++qn[2];
    else {
      qn[2] = 0;
      ++qn[0];
    }
    qn[1] = am - qn[0] - qn[2];
  }
}
#endif
#endif // 0

//
// GAMESS ordering does not yet have FOR_CART macros defined in cgshell_ordering.h, hence must handle manually
//
#if LIBINT_CGSHELL_ORDERING == LIBINT__CGSHELL__ORDERING_GAMESS
template <>
void
StdLibintTDPolicy<CGShell>::init_subobj(const StdLibintTDPolicy<CGShell>::obj_stype& cgshell,
vector<StdLibintTDPolicy<CGShell>::subobj_stype>& cgfs)
{
  const unsigned int am = TypeTraits<CGShell>::const_ref(cgshell).qn();
  unsigned int qn[3] = {0, 0, 0};

  if (am == 0) {
    subobj_stype cgf(qn);
    cgfs.push_back(cgf);
    return;
  }

  const int ammin = ((int)am + 2)/3;
  for(int am1=am; am1>=ammin; --am1) {

    for(int xyz1=0; xyz1<3; ++xyz1) {

      qn[xyz1] = am1;

      // distribute the remaining quanta according to the following rules
 
      // "nothing to distribute" is a special case
      if(am - am1 == 0) {
	std::pair<int,int> xyz(notxyz(xyz1));
	qn[xyz.first] = 0;
	qn[xyz.second] = 0;
	subobj_stype cgf(qn);
	cgfs.push_back(cgf);
      }
      else {
	int am23 = (int)am - qn[xyz1];
	const int maxam23 = std::min((int)qn[xyz1],am23);
	const int minam23 = (am23 + 1)/2;
	for(int am2=maxam23; am2>=minam23; --am2) {
	  const int xyz2min = (am2 == qn[xyz1]) ? xyz1+1 : 0;
	  for(int xyz2=xyz2min; xyz2<3; ++xyz2) {
	    if (xyz1 == xyz2)
	      continue;
	    qn[xyz2] = am2;
	    const int xyz3 = notxyz(xyz1,xyz2);
	    qn[xyz3] = am23 - am2;
	    if (qn[xyz3] == qn[xyz1] && xyz3 < xyz1 ||
		qn[xyz3] == qn[xyz2] && xyz3 < xyz2)
	      continue;
	    {
	      subobj_stype cgf(qn);
	      cgfs.push_back(cgf);
	    }
	  }
	}
      }

    }
  }

}
#else
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
    cgfs.push_back(cgf);
  END_FOR_CART
}
#endif

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
  }

  std::pair<int,int> notxyz(int a) {
    switch(a) {
    case 0: return make_pair(1,2); break;
    case 1: return make_pair(0,2); break;
    case 2: return make_pair(0,1); break;
    }
  }
}
