/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint compiler.
 *
 *  Libint compiler is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint compiler is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint compiler.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_bin_libint_gaussoper_h_
#define _libint2_src_bin_libint_gaussoper_h_

#include <braket.h>
#include <global_macros.h>
#include <prefactors.h>

#include <boost/type_traits/is_same.hpp>

namespace libint2 {

/// Applies R12vec_dot_Nabla1 to a physicists' braket
template <class F, BraketType BKType>
LinearCombination<std::shared_ptr<DGVertex>, BraketPair<F, BKType> >
R12vec_dot_Nabla1(const BraketPair<F, BKType>& bkt) {
  if (BKType == CBra || BKType == CKet)
    throw std::logic_error(
        "R12vec_dot_Nabla1 can only be applied to physicists brakets");

  const char* zeta = (BKType == PBra) ? "zeta_A" : "zeta_B";
  const char* XY = (BKType == PBra) ? "AC" : "BD";

  const auto& f = bkt[0];
  const auto& g = bkt[1];
  typedef LinearCombination<std::shared_ptr<DGVertex>, BraketPair<F, BKType> >
      ResultType;
  ResultType result;

  using namespace libint2::prefactor;
  if (f.norm())
    result += make_pair(Scalar((double)f.norm()), BraketPair<F, BKType>(f, g));

  const unsigned int nxyz = boost::is_same<F, CGF>::value ? 3 : 1;
  for (unsigned int xyz = 0; xyz < nxyz; ++xyz) {
    using namespace libint2::algebra;
    F _1 = unit<F>(xyz);
    const F& fm1 = f - _1;
    const F& gp1 = g + _1;
    if (exists(fm1)) {
      const double f_xyz = (double)(f.qn(xyz));
      result +=
          make_pair(Scalar(-1.0 * f_xyz), BraketPair<F, BKType>(fm1, gp1));
      result += make_pair(Scalar(f_xyz) * Vector(XY)[xyz],
                          BraketPair<F, BKType>(fm1, g));
    }
    const F& fp1 = f + _1;
    const F& fp2 = fp1 + _1;
    result +=
        make_pair(Scalar(-2.0) * Scalar(zeta), BraketPair<F, BKType>(fp2, g));
    result +=
        make_pair(Scalar(2.0) * Scalar(zeta), BraketPair<F, BKType>(fp1, gp1));
    result += make_pair(Scalar(-2.0) * Scalar(zeta) * Vector(XY)[xyz],
                        BraketPair<F, BKType>(fp1, g));
  }

  return result;
}

/// Applies R12vec_dot_Nabla2 to a physicists' braket
template <class F, BraketType BKType>
LinearCombination<std::shared_ptr<DGVertex>, BraketPair<F, BKType> >
R12vec_dot_Nabla2(const BraketPair<F, BKType>& bkt) {
  if (BKType == CBra || BKType == CKet)
    throw std::logic_error(
        "R12vec_dot_Nabla2 can only be applied to physicists brakets");

  const char* zeta = (BKType == PBra) ? "zeta_C" : "zeta_D";
  const char* XY = (BKType == PBra) ? "AC" : "BD";

  const F& f = bkt[0];
  const F& g = bkt[1];
  typedef LinearCombination<std::shared_ptr<DGVertex>, BraketPair<F, BKType> >
      ResultType;
  ResultType result;

  using namespace libint2::prefactor;
  if (g.norm())
    result += make_pair(Scalar(-1.0 * g.norm()), BraketPair<F, BKType>(f, g));

  const unsigned int nxyz = boost::is_same<F, CGF>::value ? 3 : 1;
  for (unsigned int xyz = 0; xyz < nxyz; ++xyz) {
    using namespace libint2::algebra;
    F _1 = unit<F>(xyz);
    const F& fp1 = f + _1;
    const F& gm1 = g - _1;
    if (exists(gm1)) {
      const double g_xyz = (double)(g.qn(xyz));
      result += make_pair(Scalar(g_xyz), BraketPair<F, BKType>(fp1, gm1));
      result += make_pair(Scalar(g_xyz) * Vector(XY)[xyz],
                          BraketPair<F, BKType>(f, gm1));
    }
    const F& gp1 = g + _1;
    const F& gp2 = gp1 + _1;
    result +=
        make_pair(Scalar(-2.0) * Scalar(zeta), BraketPair<F, BKType>(fp1, gp1));
    result +=
        make_pair(Scalar(2.0) * Scalar(zeta), BraketPair<F, BKType>(f, gp2));
    result += make_pair(Scalar(-2.0) * Scalar(zeta) * Vector(XY)[xyz],
                        BraketPair<F, BKType>(f, gp1));
  }

  return result;
}

/// Applies Nabla1 to a physicists' braket
template <class F, BraketType BKType>
LinearCombination<std::shared_ptr<DGVertex>, BraketPair<F, BKType> > Nabla1(
    const BraketPair<F, BKType>& bkt, int xyz) {
  if (BKType == CBra || BKType == CKet)
    throw std::logic_error("Nabla1 can only be applied to physicists brakets");

  const char* zeta = (BKType == PBra) ? "zeta_A" : "zeta_B";

  const F& f = bkt[0];
  const F& g = bkt[1];
  typedef LinearCombination<std::shared_ptr<DGVertex>, BraketPair<F, BKType> >
      ResultType;
  ResultType result;

  using namespace libint2::prefactor;
  using namespace libint2::algebra;
  // if applying to shell pair, change angular momentum along 0
  const unsigned int dir = boost::is_same<F, CGF>::value ? xyz : 0;
  F _1 = unit<F>(dir);
  const F& fm1 = f - _1;
  if (exists(fm1)) {
    const double f_xyz = (double)(f.qn(dir));
    result += make_pair(Scalar(f_xyz), BraketPair<F, BKType>(fm1, g));
  }
  const F& fp1 = f + _1;
  result +=
      make_pair(Scalar(-2.0) * Scalar(zeta), BraketPair<F, BKType>(fp1, g));
  return result;
}

/// Applies Nabla2 to a physicists' braket
template <class F, BraketType BKType>
LinearCombination<std::shared_ptr<DGVertex>, BraketPair<F, BKType> > Nabla2(
    const BraketPair<F, BKType>& bkt, int xyz) {
  if (BKType == CBra || BKType == CKet)
    throw std::logic_error("Nabla2 can only be applied to physicists brakets");

  const char* zeta = (BKType == PBra) ? "zeta_C" : "zeta_D";

  const F& f = bkt[0];
  const F& g = bkt[1];
  typedef LinearCombination<std::shared_ptr<DGVertex>, BraketPair<F, BKType> >
      ResultType;
  ResultType result;

  using namespace libint2::prefactor;
  using namespace libint2::algebra;
  // if applying to shell pair, change angular momentum along 0
  const unsigned int dir = boost::is_same<F, CGF>::value ? xyz : 0;
  F _1 = unit<F>(dir);
  const F& gm1 = g - _1;
  if (exists(gm1)) {
    const double g_xyz = (double)(g.qn(dir));
    result += make_pair(Scalar(g_xyz), BraketPair<F, BKType>(f, gm1));
  }
  const F& gp1 = g + _1;
  result +=
      make_pair(Scalar(-2.0) * Scalar(zeta), BraketPair<F, BKType>(f, gp1));
  return result;
}

/// Applies R12v to a physicists' braket
template <class F, BraketType BKType>
LinearCombination<std::shared_ptr<DGVertex>, BraketPair<F, BKType> > R12v(
    const BraketPair<F, BKType>& bkt, unsigned int xyz) {
  if (BKType == CBra || BKType == CKet)
    throw std::logic_error("R12v can only be applied to physicists brakets");

  const char* XY = (BKType == PBra) ? "AC" : "BD";

  const F& f = bkt[0];
  const F& g = bkt[1];
  typedef LinearCombination<std::shared_ptr<DGVertex>, BraketPair<F, BKType> >
      ResultType;
  ResultType result;

  using namespace libint2::prefactor;
  using namespace libint2::algebra;
  // if applying to shell pair, change angular momentum along 0
  const unsigned int dir = boost::is_same<F, CGF>::value ? xyz : 0;
  F _1 = unit<F>(dir);
  const F& fp1 = f + _1;
  const F& gp1 = g + _1;
  result += make_pair(Scalar(1.0), BraketPair<F, BKType>(fp1, g));
  result += make_pair(Scalar(-1.0), BraketPair<F, BKType>(f, gp1));
  result += make_pair(Vector(XY)[xyz], BraketPair<F, BKType>(f, g));
  return result;
}

};  // namespace libint2

#endif  // header guard
