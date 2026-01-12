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

#ifndef _libint2_src_bin_libint_masterintslist_h_
#define _libint2_src_bin_libint_masterintslist_h_

// need extra-long mpl list
#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_LIST_SIZE 50
#include <bfset.h>
#include <oper.h>

#include <boost/mpl/list.hpp>
#if LIBINT_SUPPORT_ONEBODYINTS
#include <integral_1_1.h>
#endif
#include <integral_11_11.h>

namespace libint2 {

//////////////////////////
// one-electron integrals
//////////////////////////
#if LIBINT_SUPPORT_ONEBODYINTS
typedef GenIntegralSet_1_1<CGShell, OverlapOper, EmptySet> Overlap_1_1_sh;
typedef GenIntegralSet_1_1<CGF, OverlapOper, EmptySet> Overlap_1_1_int;
typedef GenIntegralSet_1_1<CGShell, KineticOper, EmptySet> Kinetic_1_1_sh;
typedef GenIntegralSet_1_1<CGF, KineticOper, EmptySet> Kinetic_1_1_int;
typedef GenIntegralSet_1_1<CGShell, ElecPotOper, mType> ElecPot_1_1_sh;
typedef GenIntegralSet_1_1<CGF, ElecPotOper, mType> ElecPot_1_1_int;
typedef GenIntegralSet_1_1<CGShell, σpVσpOper, EmptySet> σpVσp_1_1_sh;
typedef GenIntegralSet_1_1<CGF, σpVσpOper, EmptySet> σpVσp_1_1_int;
typedef GenIntegralSet_1_1<CGShell, CartesianMultipoleOper<3u>, EmptySet>
    CMultipole_1_1_sh;
typedef GenIntegralSet_1_1<CGF, CartesianMultipoleOper<3u>, EmptySet>
    CMultipole_1_1_int;
typedef GenIntegralSet_1_1<CGShell, SphericalMultipoleOper, EmptySet>
    SMultipole_1_1_sh;
typedef GenIntegralSet_1_1<CGF, SphericalMultipoleOper, EmptySet>
    SMultipole_1_1_int;

typedef GenIntegralSet_1_1<CGShell1d<CartesianAxis_X>, OverlapOper, EmptySet>
    Overlap_1_1_sh_x;
typedef GenIntegralSet_1_1<CGShell1d<CartesianAxis_Y>, OverlapOper, EmptySet>
    Overlap_1_1_sh_y;
typedef GenIntegralSet_1_1<CGShell1d<CartesianAxis_Z>, OverlapOper, EmptySet>
    Overlap_1_1_sh_z;
typedef GenIntegralSet_1_1<CGF1d<CartesianAxis_X>, OverlapOper, EmptySet>
    Overlap_1_1_int_x;
typedef GenIntegralSet_1_1<CGF1d<CartesianAxis_Y>, OverlapOper, EmptySet>
    Overlap_1_1_int_y;
typedef GenIntegralSet_1_1<CGF1d<CartesianAxis_Z>, OverlapOper, EmptySet>
    Overlap_1_1_int_z;
typedef GenIntegralSet_1_1<CGShell1d<CartesianAxis_X>, KineticOper, EmptySet>
    Kinetic_1_1_sh_x;
typedef GenIntegralSet_1_1<CGShell1d<CartesianAxis_Y>, KineticOper, EmptySet>
    Kinetic_1_1_sh_y;
typedef GenIntegralSet_1_1<CGShell1d<CartesianAxis_Z>, KineticOper, EmptySet>
    Kinetic_1_1_sh_z;
typedef GenIntegralSet_1_1<CGF1d<CartesianAxis_X>, KineticOper, EmptySet>
    Kinetic_1_1_int_x;
typedef GenIntegralSet_1_1<CGF1d<CartesianAxis_Y>, KineticOper, EmptySet>
    Kinetic_1_1_int_y;
typedef GenIntegralSet_1_1<CGF1d<CartesianAxis_Z>, KineticOper, EmptySet>
    Kinetic_1_1_int_z;

typedef GenIntegralSet_1_1<CGShell1d<CartesianAxis_X>,
                           CartesianMultipoleOper<1u>, EmptySet>
    CMultipole_1_1_sh_x;
typedef GenIntegralSet_1_1<CGShell1d<CartesianAxis_Y>,
                           CartesianMultipoleOper<1u>, EmptySet>
    CMultipole_1_1_sh_y;
typedef GenIntegralSet_1_1<CGShell1d<CartesianAxis_Z>,
                           CartesianMultipoleOper<1u>, EmptySet>
    CMultipole_1_1_sh_z;
typedef GenIntegralSet_1_1<CGF1d<CartesianAxis_X>, CartesianMultipoleOper<1u>,
                           EmptySet>
    CMultipole_1_1_int_x;
typedef GenIntegralSet_1_1<CGF1d<CartesianAxis_Y>, CartesianMultipoleOper<1u>,
                           EmptySet>
    CMultipole_1_1_int_y;
typedef GenIntegralSet_1_1<CGF1d<CartesianAxis_Z>, CartesianMultipoleOper<1u>,
                           EmptySet>
    CMultipole_1_1_int_z;
#endif

//////////////////////////
// two-electron integrals
//////////////////////////
typedef GenIntegralSet_11_11<CGShell, TwoPRep, mType> TwoPRep_11_11_sq;
typedef GenIntegralSet_11_11<CGF, TwoPRep, mType> TwoPRep_11_11_int;
typedef GenIntegralSet_11_11<CGShell, R12kG12, mType> R12kG12_11_11_sq;
typedef GenIntegralSet_11_11<CGF, R12kG12, mType> R12kG12_11_11_int;
typedef GenIntegralSet_11_11<CGShell, R12kR12lG12, EmptySet>
    R12kR12lG12_11_11_sq;
typedef GenIntegralSet_11_11<CGF, R12kR12lG12, EmptySet> R12kR12lG12_11_11_int;
typedef GenIntegralSet_11_11<CGShell, TiG12, mType> TiG12_11_11_sq;
typedef GenIntegralSet_11_11<CGF, TiG12, mType> TiG12_11_11_int;
typedef GenIntegralSet_11_11<CGShell, G12TiG12, mType> G12TiG12_11_11_sq;
typedef GenIntegralSet_11_11<CGF, G12TiG12, mType> G12TiG12_11_11_int;
typedef GenIntegralSet_11_11<CGShell, DivG12prime_xTx, mType>
    DivG12prime_xTx_11_11_sq;
typedef GenIntegralSet_11_11<CGF, DivG12prime_xTx, mType>
    DivG12prime_xTx_11_11_int;
typedef GenIntegralSet<GenMultSymm2BodyOper, IncableBFSet,
                       DefaultTwoPBraket<CGShell>::Result,
                       DefaultTwoPBraket<CGShell>::Result, EmptySet>
    DummySymmIntegral_11_11_sq;
typedef GenIntegralSet<GenMultSymm2BodyOper, IncableBFSet,
                       DefaultTwoPBraket<CGF>::Result,
                       DefaultTwoPBraket<CGF>::Result, EmptySet>
    DummySymmIntegral_11_11_int;

/** All known types go into this typelist
    Every type must have a corresponding instantiation of MasterStrategy in
   strategy.cc
  */
typedef boost::mpl::list<
#if LIBINT_SUPPORT_ONEBODYINTS
    Overlap_1_1_sh, Overlap_1_1_int, Overlap_1_1_sh_x, Overlap_1_1_int_x,
    Overlap_1_1_sh_y, Overlap_1_1_int_y, Overlap_1_1_sh_z, Overlap_1_1_int_z,
    Kinetic_1_1_sh, Kinetic_1_1_int, Kinetic_1_1_sh_x, Kinetic_1_1_int_x,
    Kinetic_1_1_sh_y, Kinetic_1_1_int_y, Kinetic_1_1_sh_z, Kinetic_1_1_int_z,
    ElecPot_1_1_sh, ElecPot_1_1_int, σpVσp_1_1_sh, σpVσp_1_1_int,
    CMultipole_1_1_sh, CMultipole_1_1_int, CMultipole_1_1_sh_x,
    CMultipole_1_1_sh_y, CMultipole_1_1_sh_z, CMultipole_1_1_int_x,
    CMultipole_1_1_int_y, CMultipole_1_1_int_z, SMultipole_1_1_sh,
    SMultipole_1_1_int,
#endif
    TwoPRep_11_11_sq, TwoPRep_11_11_int, R12kG12_11_11_sq, R12kG12_11_11_int,
    R12kR12lG12_11_11_sq, R12kR12lG12_11_11_int, TiG12_11_11_sq,
    TiG12_11_11_int, G12TiG12_11_11_sq, G12TiG12_11_11_int,
    DivG12prime_xTx_11_11_sq, DivG12prime_xTx_11_11_int,
    DummySymmIntegral_11_11_sq, DummySymmIntegral_11_11_int>
    MasterIntegralTypeList;

};  // namespace libint2

#endif  // header guard
