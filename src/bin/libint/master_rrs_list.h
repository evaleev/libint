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

#ifndef _libint2_src_bin_libint_masterrrslist_h_
#define _libint2_src_bin_libint_masterrrslist_h_

#include <comp_11_DivG12prime_xTx_11.h>
#include <comp_11_g12tig12_11.h>
#include <comp_11_r12kr12lg12_11.h>
#include <comp_11_tig12_11.h>
#include <comp_1_σpVσp_1.h>
#include <comp_deriv_gauss.h>
#include <comp_xyz.h>
#include <generic_rr.h>
#include <hrr.h>
#include <itr_11_twoprep_11.h>
#include <libint2/config.h>
#include <master_ints_list.h>
#include <vrr_11_r12kg12_11.h>
#include <vrr_11_twoprep_11.h>
#include <vrr_1_onep_1.h>

#include <boost/mpl/list.hpp>

// master list of types and typedefs that describe capabilities of Libint2
// should not be used unless absolutely necessary.
// exception: strategy.cc and build_libint.cc, which depend on the entire
// typelist

namespace libint2 {

// skip derivatives with respect to B in (A|B)
constexpr auto trinvskip1_part = 0;
constexpr auto trinvskip1_where = InKet;
// skip derivatives with respect to C or D in (AB|CD)
constexpr auto trinvskip2_part = 1;
constexpr auto trinvskip2_where =
    (LIBINT_SHELL_SET == LIBINT_SHELL_SET_ORCA) ? InKet : InBra;

///////////
// RRs
///////////
typedef HRR<TwoPRep_11_11_sq, CGShell, 0, InBra, 0, InKet, 0>
    HRR_ab_11_TwoPRep_11_sh;
typedef HRR<TwoPRep_11_11_sq, CGShell, 1, InBra, 0, InKet, 0>
    HRR_cd_11_TwoPRep_11_sh;
typedef HRR<TwoPRep_11_11_sq, CGShell, 0, InKet, 0, InBra, 0>
    HRR_ba_11_TwoPRep_11_sh;
typedef HRR<TwoPRep_11_11_sq, CGShell, 1, InKet, 0, InBra, 0>
    HRR_dc_11_TwoPRep_11_sh;
typedef HRR<R12kG12_11_11_sq, CGShell, 0, InBra, 0, InKet, 0>
    HRR_ab_11_R12kG12_11_sh;
typedef HRR<R12kG12_11_11_sq, CGShell, 1, InBra, 0, InKet, 0>
    HRR_cd_11_R12kG12_11_sh;
typedef HRR<R12kG12_11_11_sq, CGShell, 0, InKet, 0, InBra, 0>
    HRR_ba_11_R12kG12_11_sh;
typedef HRR<R12kG12_11_11_sq, CGShell, 1, InKet, 0, InBra, 0>
    HRR_dc_11_R12kG12_11_sh;
typedef HRR<R12kR12lG12_11_11_sq, CGShell, 0, InBra, 0, InKet, 0>
    HRR_ab_11_R12kR12lG12_11_sh;
typedef HRR<R12kR12lG12_11_11_sq, CGShell, 1, InBra, 0, InKet, 0>
    HRR_cd_11_R12kR12lG12_11_sh;
typedef HRR<TiG12_11_11_sq, CGShell, 0, InBra, 0, InKet, 0>
    HRR_ab_11_TiG12_11_sh;
typedef HRR<TiG12_11_11_sq, CGShell, 1, InBra, 0, InKet, 0>
    HRR_cd_11_TiG12_11_sh;
typedef HRR<G12TiG12_11_11_sq, CGShell, 0, InBra, 0, InKet, 0>
    HRR_ab_11_G12TiG12_11_sh;
typedef HRR<G12TiG12_11_11_sq, CGShell, 1, InBra, 0, InKet, 0>
    HRR_cd_11_G12TiG12_11_sh;
typedef HRR<DivG12prime_xTx_11_11_sq, CGShell, 0, InBra, 0, InKet, 0>
    HRR_ab_11_DivG12prime_xTx_sh;
typedef HRR<DivG12prime_xTx_11_11_sq, CGShell, 1, InBra, 0, InKet, 0>
    HRR_cd_11_DivG12prime_xTx_sh;

typedef HRR<TwoPRep_11_11_int, CGF, 0, InBra, 0, InKet, 0>
    HRR_ab_11_TwoPRep_11_int;
typedef HRR<TwoPRep_11_11_int, CGF, 1, InBra, 0, InKet, 0>
    HRR_cd_11_TwoPRep_11_int;
typedef HRR<TwoPRep_11_11_int, CGF, 0, InKet, 0, InBra, 0>
    HRR_ba_11_TwoPRep_11_int;
typedef HRR<TwoPRep_11_11_int, CGF, 1, InKet, 0, InBra, 0>
    HRR_dc_11_TwoPRep_11_int;
typedef HRR<R12kG12_11_11_int, CGF, 0, InBra, 0, InKet, 0>
    HRR_ab_11_R12kG12_11_int;
typedef HRR<R12kG12_11_11_int, CGF, 1, InBra, 0, InKet, 0>
    HRR_cd_11_R12kG12_11_int;
typedef HRR<R12kG12_11_11_int, CGF, 0, InKet, 0, InBra, 0>
    HRR_ba_11_R12kG12_11_int;
typedef HRR<R12kG12_11_11_int, CGF, 1, InKet, 0, InBra, 0>
    HRR_dc_11_R12kG12_11_int;
typedef HRR<R12kR12lG12_11_11_int, CGF, 0, InBra, 0, InKet, 0>
    HRR_ab_11_R12kR12lG12_11_int;
typedef HRR<R12kR12lG12_11_11_int, CGF, 1, InBra, 0, InKet, 0>
    HRR_cd_11_R12kR12lG12_11_int;
typedef HRR<TiG12_11_11_int, CGF, 0, InBra, 0, InKet, 0> HRR_ab_11_TiG12_11_int;
typedef HRR<TiG12_11_11_int, CGF, 1, InBra, 0, InKet, 0> HRR_cd_11_TiG12_11_int;
typedef HRR<G12TiG12_11_11_int, CGF, 0, InBra, 0, InKet, 0>
    HRR_ab_11_G12TiG12_11_int;
typedef HRR<G12TiG12_11_11_int, CGF, 1, InBra, 0, InKet, 0>
    HRR_cd_11_G12TiG12_11_int;
typedef HRR<DivG12prime_xTx_11_11_int, CGF, 0, InBra, 0, InKet, 0>
    HRR_ab_11_DivG12prime_xTx_int;
typedef HRR<DivG12prime_xTx_11_11_int, CGF, 1, InBra, 0, InKet, 0>
    HRR_cd_11_DivG12prime_xTx_int;

typedef HRR<DummySymmIntegral_11_11_sq, CGShell, 0, InBra, 0, InKet, 0>
    HRR_ab_11_Dummy_11_sh;
typedef HRR<DummySymmIntegral_11_11_sq, CGShell, 1, InBra, 0, InKet, 0>
    HRR_cd_11_Dummy_11_sh;
typedef HRR<DummySymmIntegral_11_11_sq, CGShell, 0, InKet, 0, InBra, 0>
    HRR_ba_11_Dummy_11_sh;
typedef HRR<DummySymmIntegral_11_11_sq, CGShell, 1, InKet, 0, InBra, 0>
    HRR_dc_11_Dummy_11_sh;
typedef HRR<DummySymmIntegral_11_11_int, CGF, 0, InBra, 0, InKet, 0>
    HRR_ab_11_Dummy_11_int;
typedef HRR<DummySymmIntegral_11_11_int, CGF, 1, InBra, 0, InKet, 0>
    HRR_cd_11_Dummy_11_int;
typedef HRR<DummySymmIntegral_11_11_int, CGF, 0, InKet, 0, InBra, 0>
    HRR_ba_11_Dummy_11_int;
typedef HRR<DummySymmIntegral_11_11_int, CGF, 1, InKet, 0, InBra, 0>
    HRR_dc_11_Dummy_11_int;

#if LIBINT_SUPPORT_ONEBODYINTS
typedef HRR<Overlap_1_1_sh, CGShell, 0, InBra, 0, InKet, 0>
    HRR_ab_1_Overlap_1_sh;
typedef HRR<Overlap_1_1_int, CGF, 0, InBra, 0, InKet, 0> HRR_ab_1_Overlap_1_int;
typedef HRR<Overlap_1_1_sh, CGShell, 0, InKet, 0, InBra, 0>
    HRR_ba_1_Overlap_1_sh;
typedef HRR<Overlap_1_1_int, CGF, 0, InKet, 0, InBra, 0> HRR_ba_1_Overlap_1_int;

typedef HRR<ElecPot_1_1_sh, CGShell, 0, InBra, 0, InKet, 0>
    HRR_ab_1_ElecPot_1_sh;
typedef HRR<ElecPot_1_1_int, CGF, 0, InBra, 0, InKet, 0> HRR_ab_1_ElecPot_1_int;
typedef HRR<ElecPot_1_1_sh, CGShell, 0, InKet, 0, InBra, 0>
    HRR_ba_1_ElecPot_1_sh;
typedef HRR<ElecPot_1_1_int, CGF, 0, InKet, 0, InBra, 0> HRR_ba_1_ElecPot_1_int;

typedef VRR_1_Overlap_1<CGShell, InBra> VRR_a_1_Overlap_1_sh;
typedef VRR_1_Overlap_1<CGF, InBra> VRR_a_1_Overlap_1_int;
typedef VRR_1_Overlap_1<CGShell, InKet> VRR_b_1_Overlap_1_sh;
typedef VRR_1_Overlap_1<CGF, InKet> VRR_b_1_Overlap_1_int;

typedef VRR_1_Overlap_1_1d<CartesianAxis_X, InBra> VRR_a_1_Overlap_1_int_x;
typedef VRR_1_Overlap_1_1d<CartesianAxis_Y, InBra> VRR_a_1_Overlap_1_int_y;
typedef VRR_1_Overlap_1_1d<CartesianAxis_Z, InBra> VRR_a_1_Overlap_1_int_z;
typedef VRR_1_Overlap_1_1d<CartesianAxis_X, InKet> VRR_b_1_Overlap_1_int_x;
typedef VRR_1_Overlap_1_1d<CartesianAxis_Y, InKet> VRR_b_1_Overlap_1_int_y;
typedef VRR_1_Overlap_1_1d<CartesianAxis_Z, InKet> VRR_b_1_Overlap_1_int_z;

typedef VRR_1_Kinetic_1<CGShell, InBra> VRR_a_1_Kinetic_1_sh;
typedef VRR_1_Kinetic_1<CGF, InBra> VRR_a_1_Kinetic_1_int;
typedef VRR_1_Kinetic_1<CGShell, InKet> VRR_b_1_Kinetic_1_sh;
typedef VRR_1_Kinetic_1<CGF, InKet> VRR_b_1_Kinetic_1_int;

typedef VRR_1_ElecPot_1<CGShell, InBra> VRR_a_1_ElecPot_1_sh;
typedef VRR_1_ElecPot_1<CGF, InBra> VRR_a_1_ElecPot_1_int;
typedef VRR_1_ElecPot_1<CGShell, InKet> VRR_b_1_ElecPot_1_sh;
typedef VRR_1_ElecPot_1<CGF, InKet> VRR_b_1_ElecPot_1_int;

typedef CR_1_σpVσp_1<CGShell> CR_1_σpVσp_1_sh;
typedef CR_1_σpVσp_1<CGF> CR_1_σpVσp_1_int;

// TODO investigate whether need to stay away from HRR for now to be sure that
// multipoles are computed as precisely as possible
typedef HRR<SMultipole_1_1_sh, CGShell, 0, InBra, 0, InKet, 0>
    HRR_ab_1_SMultipole_1_sh;
typedef HRR<SMultipole_1_1_int, CGF, 0, InBra, 0, InKet, 0>
    HRR_ab_1_SMultipole_1_int;
typedef HRR<SMultipole_1_1_sh, CGShell, 0, InKet, 0, InBra, 0>
    HRR_ba_1_SMultipole_1_sh;
typedef HRR<SMultipole_1_1_int, CGF, 0, InKet, 0, InBra, 0>
    HRR_ba_1_SMultipole_1_int;
typedef VRR_1_SMultipole_1<CGShell, InBra> VRR_a_1_SMultipole_1_sh;
typedef VRR_1_SMultipole_1<CGF, InBra> VRR_a_1_SMultipole_1_int;
typedef VRR_1_SMultipole_1<CGShell, InKet> VRR_b_1_SMultipole_1_sh;
typedef VRR_1_SMultipole_1<CGF, InKet> VRR_b_1_SMultipole_1_int;
#endif

typedef VRR_11_TwoPRep_11<CGShell, 0, InBra> VRR_a_11_TwoPRep_11_sh;
typedef VRR_11_TwoPRep_11<CGShell, 1, InBra> VRR_c_11_TwoPRep_11_sh;
typedef VRR_11_TwoPRep_11<CGF, 0, InBra> VRR_a_11_TwoPRep_11_int;
typedef VRR_11_TwoPRep_11<CGF, 1, InBra> VRR_c_11_TwoPRep_11_int;
typedef VRR_11_TwoPRep_11<CGShell, 0, InKet> VRR_b_11_TwoPRep_11_sh;
typedef VRR_11_TwoPRep_11<CGShell, 1, InKet> VRR_d_11_TwoPRep_11_sh;
typedef VRR_11_TwoPRep_11<CGF, 0, InKet> VRR_b_11_TwoPRep_11_int;
typedef VRR_11_TwoPRep_11<CGF, 1, InKet> VRR_d_11_TwoPRep_11_int;

typedef VRR_11_R12kG12_11<CGShell, 0, InBra> VRR_a_11_R12kG12_11_sh;
typedef VRR_11_R12kG12_11<CGShell, 1, InBra> VRR_c_11_R12kG12_11_sh;
typedef VRR_11_R12kG12_11<CGF, 0, InBra> VRR_a_11_R12kG12_11_int;
typedef VRR_11_R12kG12_11<CGF, 1, InBra> VRR_c_11_R12kG12_11_int;
typedef VRR_11_R12kG12_11<CGShell, 0, InKet> VRR_b_11_R12kG12_11_sh;
typedef VRR_11_R12kG12_11<CGShell, 1, InKet> VRR_d_11_R12kG12_11_sh;
typedef VRR_11_R12kG12_11<CGF, 0, InKet> VRR_b_11_R12kG12_11_int;
typedef VRR_11_R12kG12_11<CGF, 1, InKet> VRR_d_11_R12kG12_11_int;

typedef CR_11_R12kR12lG12_11<CGShell> CR_11_R12kR12lG12_11_sh;
typedef CR_11_R12kR12lG12_11<CGF> CR_11_R12kR12lG12_11_int;

typedef CR_11_TiG12_11<CGShell> CR_11_TiG12_11_sh;
typedef CR_11_TiG12_11<CGF> CR_11_TiG12_11_int;

typedef CR_11_G12TiG12_11<CGShell> CR_11_G12TiG12_11_sh;
typedef CR_11_G12TiG12_11<CGF> CR_11_G12TiG12_11_int;

typedef CR_11_DivG12prime_xTx_11<CGShell> CR_11_DivG12prime_xTx_11_sh;
typedef CR_11_DivG12prime_xTx_11<CGF> CR_11_DivG12prime_xTx_11_int;

typedef ITR_11_TwoPRep_11<GenIntegralSet_11_11, CGShell, 0, InBra>
    ITR_a_11_TwoPRep_11_sh;
typedef ITR_11_TwoPRep_11<GenIntegralSet_11_11, CGF, 0, InBra>
    ITR_a_11_TwoPRep_11_int;
typedef ITR_11_TwoPRep_11<GenIntegralSet_11_11, CGShell, 0, InKet>
    ITR_b_11_TwoPRep_11_sh;
typedef ITR_11_TwoPRep_11<GenIntegralSet_11_11, CGF, 0, InKet>
    ITR_b_11_TwoPRep_11_int;
typedef ITR_11_TwoPRep_11<GenIntegralSet_11_11, CGShell, 1, InBra>
    ITR_c_11_TwoPRep_11_sh;
typedef ITR_11_TwoPRep_11<GenIntegralSet_11_11, CGF, 1, InBra>
    ITR_c_11_TwoPRep_11_int;
typedef ITR_11_TwoPRep_11<GenIntegralSet_11_11, CGShell, 1, InKet>
    ITR_d_11_TwoPRep_11_sh;
typedef ITR_11_TwoPRep_11<GenIntegralSet_11_11, CGF, 1, InKet>
    ITR_d_11_TwoPRep_11_int;

typedef CR_DerivGauss<TwoPRep_11_11_sq, 0, InBra, trinvskip2_part,
                      trinvskip2_where>
    Deriv_a_11_TwoPRep_11_sh;
typedef CR_DerivGauss<TwoPRep_11_11_sq, 0, InKet, trinvskip2_part,
                      trinvskip2_where>
    Deriv_b_11_TwoPRep_11_sh;
typedef CR_DerivGauss<TwoPRep_11_11_sq, 1, InBra, trinvskip2_part,
                      trinvskip2_where>
    Deriv_c_11_TwoPRep_11_sh;
typedef CR_DerivGauss<TwoPRep_11_11_sq, 1, InKet, trinvskip2_part,
                      trinvskip2_where>
    Deriv_d_11_TwoPRep_11_sh;
typedef CR_DerivGauss<TwoPRep_11_11_int, 0, InBra, trinvskip2_part,
                      trinvskip2_where>
    Deriv_a_11_TwoPRep_11_int;
typedef CR_DerivGauss<TwoPRep_11_11_int, 0, InKet, trinvskip2_part,
                      trinvskip2_where>
    Deriv_b_11_TwoPRep_11_int;
typedef CR_DerivGauss<TwoPRep_11_11_int, 1, InBra, trinvskip2_part,
                      trinvskip2_where>
    Deriv_c_11_TwoPRep_11_int;
typedef CR_DerivGauss<TwoPRep_11_11_int, 1, InKet, trinvskip2_part,
                      trinvskip2_where>
    Deriv_d_11_TwoPRep_11_int;

};  // namespace libint2

#endif  // header guard
