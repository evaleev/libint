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

#ifndef _libint2_include_libint2_util_mathconstants_h_
#define _libint2_include_libint2_util_mathconstants_h_

/// @file math_constants.h
/// mathematical constants used by the C++ API headers; these replace the
/// non-standard M_PI etc. macros, which MSVC only provides if
/// _USE_MATH_DEFINES is defined before <cmath> is included

namespace libint2 {
namespace math {

/// \f$ \pi \f$
constexpr double pi = 3.141592653589793238462643383279502884;
/// \f$ \sqrt{2} \f$
constexpr double sqrt2 = 1.414213562373095048801688724209698079;

}  // namespace math
}  // namespace libint2

#endif  // _libint2_include_libint2_util_mathconstants_h_
