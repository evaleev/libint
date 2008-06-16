
#ifndef _libint2_src_bin_libint_gaussoper_h_
#define _libint2_src_bin_libint_gaussoper_h_

#include <braket.h>
#include <prefactors.h>
#include <global_macros.h>
#include <boost/mpl/if.hpp>

namespace libint2 {

  using namespace boost;
  using namespace boost::mpl;

  /// Applies R12vec_dot_Nabla1 to a physicists' braket
  template <class F, BraketType BKType>
  LinearCombination<
    SafePtr<DGVertex>,
    BraketPair<F,BKType>
  > R12vec_dot_Nabla1(const BraketPair<F,BKType>& bkt) {
    if (BKType == CBra || BKType == CKet)
      throw std::logic_error("R12vec_dot_Nabla can only be applied to physicists brakets");

    const char* twozeta = (BKType == CBra) ? "twozeta_a" : "twozeta_b";
    const char* XY = (BKType == CBra) ? "AC" : "BD";
    
    const F& f = bkt[0];
    const F& g = bkt[1];
    typedef LinearCombination< SafePtr<DGVertex> ,BraketPair<F,BKType> > ResultType;
    ResultType result;
    
    using namespace libint2::prefactor;
    if (f.norm())
      result += make_pair(Scalar(f.norm()),
                          BraketPair<F,BKType>(f,g));
      
    const unsigned int nxyz = is_same<F,CGF>::value ? 3 : 1;
    for(unsigned int xyz=0; xyz<nxyz; ++xyz) {
      using namespace libint2::algebra;
      F _1 = unit<F>(xyz);
      const F& fm1 = f - _1;
      const F& gp1 = g + _1;
      if (exists(fm1)) {
        const int f_xyz = static_cast<int>(f.qn(xyz));
        result += make_pair(Scalar(-1*f_xyz),
                            BraketPair<F,BKType>(fm1,gp1));
        result += make_pair(Scalar(f_xyz)*Vector(XY)[xyz],
                            BraketPair<F,BKType>(fm1,g));
      }
      const F& fp1 = f + _1;
      const F& fp2 = fp1 + _1;
      result += make_pair(Scalar(-1) * Scalar(twozeta),
                          BraketPair<F,BKType>(fp2,g));
      result += make_pair(Scalar(twozeta),
                          BraketPair<F,BKType>(fp1,gp1));
      result += make_pair(Scalar(-1) * Scalar(twozeta) * Vector(XY)[xyz],
                          BraketPair<F,BKType>(fp1,g));
    }
    
    return result;
  }

  /// Applies R12vec_dot_Nabla2 to a physicists' braket
  template <class F, BraketType BKType>
  LinearCombination<
    SafePtr<DGVertex>,
    BraketPair<F,BKType>
  > R12vec_dot_Nabla2(const BraketPair<F,BKType>& bkt) {
    if (BKType == CBra || BKType == CKet)
      throw std::logic_error("R12vec_dot_Nabla can only be applied to physicists brakets");

    const char* twozeta = (BKType == CBra) ? "twozeta_c" : "twozeta_d";
    const char* XY = (BKType == CBra) ? "AC" : "BD";
    
    const F& f = bkt[0];
    const F& g = bkt[1];
    typedef LinearCombination< SafePtr<DGVertex> ,BraketPair<F,BKType> > ResultType;
    ResultType result;
    
    using namespace libint2::prefactor;
    if (g.norm())
      result += make_pair(Scalar(-(int)g.norm()),
                          BraketPair<F,BKType>(f,g));
      
    const unsigned int nxyz = is_same<F,CGF>::value ? 3 : 1;
    for(unsigned int xyz=0; xyz<nxyz; ++xyz) {
      using namespace libint2::algebra;
      F _1 = unit<F>(xyz);
      const F& fp1 = f + _1;
      const F& gm1 = g - _1;
      if (exists(gm1)) {
        const int g_xyz = static_cast<int>(g.qn(xyz));
        result += make_pair(Scalar(g_xyz),
                            BraketPair<F,BKType>(fp1,gm1));
        result += make_pair(Scalar(g_xyz)*Vector(XY)[xyz],
                            BraketPair<F,BKType>(f,gm1));
      }
      const F& gp1 = g + _1;
      const F& gp2 = gp1 + _1;
      result += make_pair(Scalar(twozeta),
                          BraketPair<F,BKType>(f,gp2));
      result += make_pair(Scalar(-1) * Scalar(twozeta),
                          BraketPair<F,BKType>(fp1,gp1));
      result += make_pair(Scalar(-1) * Scalar(twozeta) * Vector(XY)[xyz],
                          BraketPair<F,BKType>(f,gp1));
    }
    
    return result;
  }

};

#endif // header guard
