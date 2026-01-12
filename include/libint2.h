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

#ifndef _libint2_header_
#define _libint2_header_

#define LIBINT_T_SS_EREP_SS(mValue) \
  _aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_##mValue
#define LIBINT_T_SS_Km1G12_SS(mValue) \
  _aB_s___0__s___1___r12_minus_1_g12_s___0__s___1___Ab__up_##mValue
#define LIBINT_T_SS_K0G12_SS_0 \
  _aB_s___0__s___1___r12_0_g12_s___0__s___1___Ab__up_0
#define LIBINT_T_SS_K2G12_SS_0 \
  _aB_s___0__s___1___r12_2_g12_s___0__s___1___Ab__up_0
#define LIBINT_T_SS_K4G12_SS_0 \
  _aB_s___0__s___1___r12_4_g12_s___0__s___1___Ab__up_0
#define LIBINT_T_S_OVERLAP_S _aB_s___0___Overlap_s___0___Ab__up_
#define LIBINT_T_S_KINETIC_S _aB_s___0___Kinetic_s___0___Ab__up_
#define LIBINT_T_S_ELECPOT_S(mValue) _aB_s___0___ElecPot_s___0___Ab__up_##mValue

#include <libint2/util/configuration.h>
#include <libint2/util/generated/libint2_params.h>
#include <libint2/util/generated/libint2_types.h>
#include <libint2/util/intrinsic_types.h>

#if defined(__cplusplus)
#include <libint2/numeric.h>
#endif

#endif /* header guard */

#include <libint2/util/generated/libint2_iface.h>
