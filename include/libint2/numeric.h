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

#ifndef _libint2_include_numeric_h_
#define _libint2_include_numeric_h_

#include <libint2/config.h>

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <type_traits>

#if LIBINT_HAS_MPFR
#include <gmpxx.h>
#include <mpfr.h>

#include <cstddef>
#endif

#include <libint2/util/generated/libint2_params.h>
#include <libint2/util/type_traits.h>

#if LIBINT_HAS_MPFR
/// implement exp for mpf_class using MPFR ... I do not claim to know what
/// issues the rounding presents here
inline mpf_class exp(mpf_class x) {
  const auto prec = x.get_prec();
  mpfr_t x_r;
  mpfr_init2(x_r, prec);
  mpfr_set_f(x_r, x.get_mpf_t(), MPFR_RNDN);

  mpfr_t expx_r;
  mpfr_init2(expx_r, prec);
  mpfr_exp(expx_r, x_r, MPFR_RNDN);

  mpf_t expx;
  mpf_init2(expx, prec);
  mpfr_get_f(expx, expx_r, MPFR_RNDN);
  mpf_class result(expx, prec);

  mpfr_clear(x_r);
  mpfr_clear(expx_r);
  mpf_clear(expx);

  return result;
}
/// implement pow for mpf_class using MPFR ... I do not claim to know what
/// issues the rounding presents here
inline mpf_class pow(mpf_class x, int a) {
  const auto prec = x.get_prec();
  mpf_t x_to_a;
  mpf_init2(x_to_a, prec);
  if (a >= 0)
    mpf_pow_ui(x_to_a, x.get_mpf_t(), (unsigned int)a);
  else
    mpf_pow_ui(x_to_a, x.get_mpf_t(), (unsigned int)(-a));
  mpf_class result(x_to_a, prec);
  if (a < 0) result = 1.0 / result;
  mpf_clear(x_to_a);
  return result;
}
#ifndef _MSC_VER
/// this is needed to avoid ambiguity in pow(2.0, 2) ... the above pow competes
/// with standard double pow(double, double)
inline double pow(double x, int a) {
  return std::pow(x, static_cast<double>(a));
}
#endif
/// implement erf for mpf_class using MPFR ... I do not claim to know what
/// issues the rounding presents here
inline mpf_class erf(mpf_class x) {
  const auto prec = x.get_prec();
  mpfr_t x_r;
  mpfr_init2(x_r, prec);
  mpfr_set_f(x_r, x.get_mpf_t(), MPFR_RNDN);

  mpfr_t erfx_r;
  mpfr_init2(erfx_r, prec);
  mpfr_erf(erfx_r, x_r, MPFR_RNDN);

  mpf_t erfx;
  mpf_init2(erfx, prec);
  mpfr_get_f(erfx, erfx_r, MPFR_RNDN);
  mpf_class result(erfx, prec);

  mpfr_clear(x_r);
  mpfr_clear(erfx_r);
  mpf_clear(erfx);

  return result;
}
/// implement acos for mpf_class using MPFR ... I do not claim to know what
/// issues the rounding presents here
inline mpf_class acos(mpf_class x) {
  const auto prec = x.get_prec();
  mpfr_t x_r;
  mpfr_init2(x_r, prec);
  mpfr_set_f(x_r, x.get_mpf_t(), MPFR_RNDN);

  mpfr_t acosx_r;
  mpfr_init2(acosx_r, prec);
  mpfr_acos(acosx_r, x_r, MPFR_RNDN);

  mpf_t acosx;
  mpf_init2(acosx, prec);
  mpfr_get_f(acosx, acosx_r, MPFR_RNDN);
  mpf_class result(acosx, prec);

  mpfr_clear(x_r);
  mpfr_clear(acosx_r);
  mpf_clear(acosx);

  return result;
}
/// implement log for mpf_class using MPFR ... I do not claim to know what
/// issues the rounding presents here
inline mpf_class log(mpf_class x) {
  const auto prec = x.get_prec();
  mpfr_t x_r;
  mpfr_init2(x_r, prec);
  mpfr_set_f(x_r, x.get_mpf_t(), MPFR_RNDN);

  mpfr_t logx_r;
  mpfr_init2(logx_r, prec);
  mpfr_log(logx_r, x_r, MPFR_RNDN);

  mpf_t logx;
  mpf_init2(logx, prec);
  mpfr_get_f(logx, logx_r, MPFR_RNDN);
  mpf_class result(logx, prec);

  mpfr_clear(x_r);
  mpfr_clear(logx_r);
  mpf_clear(logx);

  return result;
}
#endif

#ifdef LIBINT_HAS_MPFR
using LIBINT2_REF_REALTYPE = mpf_class;
#else
using LIBINT2_REF_REALTYPE = double;
#endif

namespace libint2 {
using value_type = LIBINT2_REALTYPE;
using scalar_type = libint2::vector_traits<value_type>::scalar_type;

template <typename Real>
inline Real get_epsilon(const Real& value);

#ifdef LIBINT_HAS_MPFR
template <>
inline mpf_class get_epsilon(const mpf_class& value) {
  const auto nbits = value.get_prec();
  return pow(mpf_class(2, nbits), -nbits);
};
#endif

template <typename Real>
inline Real get_epsilon(const Real& value) {
  return std::numeric_limits<Real>::epsilon();
}

template <typename Real>
inline int get_max_digits10(const Real& value);

#ifdef LIBINT_HAS_MPFR
template <>
inline int get_max_digits10(const mpf_class& value) {
  const auto nbits = value.get_prec();
  return std::ceil(nbits * std::log10(2) + 1);
};
#endif

template <typename Real>
inline int get_max_digits10(const Real& value) {
  return std::numeric_limits<Real>::max_digits10;
}

template <typename To, typename From>
typename std::enable_if<!std::is_same<typename std::decay<To>::type,
                                      typename std::decay<From>::type>::value,
                        To>::type
sstream_convert(From&& from) {
  std::stringstream ss;
  ss << std::scientific << std::setprecision(get_max_digits10(from)) << from;
  To to(ss.str().c_str());
  return to;
}

template <typename To>
To sstream_convert(const To& from) {
  return from;
}

};  // namespace libint2

#endif  // _libint2_include_numeric_h_
