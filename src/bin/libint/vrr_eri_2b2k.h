
#ifndef _libint2_src_bin_libint_vrreri2b2k_h_
#define _libint2_src_bin_libint_vrreri2b2k_h_

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <rr.h>

using namespace std;


namespace libint2 {

  template <class BFSet>
    VRR_ERI_2b2k<BFSet>::VRR_ERI_2b2k(const TwoERep_2b2k<BFSet>* Tint)
    {
      target_ = Tint;

      BFSet sh_a = Tint->bra(0,0);
      BFSet sh_b = Tint->ket(0,0);
      BFSet sh_c = Tint->bra(1,0);
      BFSet sh_d = Tint->ket(1,0);
      unsigned int m = Tint->m();

      vector<BFSet> bra[2];
      vector<BFSet> ket[2];

      bra[0].push_back(sh_a);
      bra[1].push_back(sh_c);
      ket[0].push_back(sh_b);
      ket[1].push_back(sh_d);

      bra[0][0].dec();
      children_[0] = new TwoERep_2b2k<BFSet>(bra,ket,m);
      children_[1] = new TwoERep_2b2k<BFSet>(bra,ket,m+1);
      bra[0][0].dec();
      children_[2] = new TwoERep_2b2k<BFSet>(bra,ket,m);
      children_[3] = new TwoERep_2b2k<BFSet>(bra,ket,m+1);
      bra[0][0].inc();

      bra[1][0].dec();
      children_[4] = new TwoERep_2b2k<BFSet>(bra,ket,m);

    };

  template <class BFSet>
    VRR_ERI_2b2k<BFSet>::~VRR_ERI_2b2k()
    {
    };

};

#endif
