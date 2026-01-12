/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint library.
 *
 *  Libint library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_lib_libint_boysfwd_h_
#define _libint2_src_lib_libint_boysfwd_h_

namespace libint2 {

template <typename Real>
struct FmEval_Reference;

template <typename Real>
struct FmEval_Reference2;

template <typename Real>
class FmEval_Chebyshev7;

template <typename Real, int INTERPOLATION_ORDER>
class FmEval_Taylor;

template <typename Real>
struct TennoGmEval;

template <typename Real, int K>
struct GaussianGmEval;

namespace detail {
/// some evaluators need thread-local scratch, but most don't
template <typename CoreEval>
struct CoreEvalScratch;
}  // namespace detail

template <typename GmEvalFunction>
struct GenericGmEval;

/// Obara-Saika core ints code
namespace os_core_ints {
template <typename Real>
struct delta_gm_eval;
template <typename Real, int K>
struct r12_xx_K_gm_eval;
template <typename Real>
struct erf_coulomb_gm_eval;
template <typename Real>
struct erfx_coulomb_gm_eval;
}  // namespace os_core_ints

/*
 *  Slater geminal fitting is available only if have LAPACK
 */
#if HAVE_LAPACK
template <typename Real>
void stg_ng_fit(unsigned int n, Real zeta,
                std::vector<std::pair<Real, Real> >& geminal, Real xmin = 0.0,
                Real xmax = 10.0, unsigned int npts = 1001);
#endif

}  // namespace libint2

#endif  // header guard
