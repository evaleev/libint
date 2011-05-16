
#ifndef _libint2_src_bin_libint_integraldecl_h_
#define _libint2_src_bin_libint_integraldecl_h_

namespace libint2 {

  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
  class GenIntegralSet;
  template <class Oper, class BFS, class AuxQuanta>
  class GenIntegralSet_11_11;

#if 0
  template <class BFS> class TwoPRep_11_11;
#endif
  template <class BFS, int K> class R12kG12_11_11;
  template <class BFS, int K> class TiG12_11_11;
  template <class BFS> class R1dotR1G12_11_11;
  template <class BFS> class R2dotR2G12_11_11;
  template <class BFS> class R1dotR2G12_11_11;
  
};

#endif
