
#ifndef _libint2_src_bin_libint_policy_h_
#define _libint2_src_bin_libint_policy_h_

#include <rr.h>

using namespace std;


namespace libint2 {

  /** Policy template used to construct a policy which specifies all information
  about assumptions of the generated code. Among such are orderings of
  Gaussians within shells, ordering of operators within sets, etc.

  The only parameter so far, CGShellOrder, is a pointer to the function
  which generates CGFs in the canonical order for a given CGShell.

  The default policy is StdLibintPolicy.
  */
  template <void (*CGShellOrder)(const CGShell*, std::vector< const CGF* >&) >
  struct Policy {
    static void cgshell_to_cgfvector(const CGShell* cgshell, std::vector< const CGF* >& cgfs) {
      CGShellOrder(cgshell, cgfs);
    }
  };

  /** StdLibintPolicy describes assumptions about orderings, etc. in Libint version 1.

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
  struct StdLibint1Order {
    static void cgshell_to_cgfs(const CGShell* cgshell, vector<const CGF*>& cgfs)
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
  };

  typedef Policy<&StdLibint1Order::cgshell_to_cgfs> StdLibintPolicy;

};

#endif
