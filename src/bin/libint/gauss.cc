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

#include <bfset.h>
#include <ctype.h>
#include <default_params.h>
#include <exception.h>
#include <string.h>

#include <cassert>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace std;
using namespace libint2;

namespace libint2 {
LIBINT2_UINT_LEAST64 cgshell_key_l_offset_[] = {
    0,    1,    4,    10,   20,   35,   56,   84,   120,  165,  220,
    286,  364,  455,  560,  680,  816,  969,  1140, 1330, 1540, 1771,
    2024, 2300, 2600, 2925, 3276, 3654, 4060, 4495, 4960, 5456, 5984};
LIBINT2_UINT_LEAST64 oderiv_1d_key_l_offset_[] = {0, 1, 2, 3, 4};
LIBINT2_UINT_LEAST64 oderiv_2d_key_l_offset_[] = {0, 1, 3, 6, 10};
LIBINT2_UINT_LEAST64 oderiv_3d_key_l_offset_[] = {0, 1, 4, 10, 20};

template <typename T, std::size_t N>
std::array<T, N> make_std_array(T* data) {
  std::array<T, N> result;
  std::copy(data, data + N, result.begin());
  return result;
}
};  // namespace libint2

std::array<LIBINT2_UINT_LEAST64, CGShell::max_qn + 1> CGF::key_l_offset(
    make_std_array<LIBINT2_UINT_LEAST64, CGShell::max_qn + 1>(
        cgshell_key_l_offset_));

template <>
std::array<LIBINT2_UINT_LEAST64, OriginDerivative<1u>::max_deriv + 1>
    OriginDerivative<1u>::key_l_offset(
        make_std_array<LIBINT2_UINT_LEAST64, OriginDerivative::max_deriv + 1>(
            oderiv_1d_key_l_offset_));
template <>
std::array<LIBINT2_UINT_LEAST64, OriginDerivative<3u>::max_deriv + 1>
    OriginDerivative<3u>::key_l_offset(
        make_std_array<LIBINT2_UINT_LEAST64, OriginDerivative::max_deriv + 1>(
            oderiv_3d_key_l_offset_));

namespace {
std::string am_to_symbol(unsigned int l, bool contracted) {
  std::string result;
  const size_t lmax_plus_1 = sizeof(LIBINT_AM2SYMBOL) - 1;
  do {
    const unsigned int digit = l % lmax_plus_1;
    char letter = LIBINT_AM2SYMBOL[digit];
    if (contracted) letter = toupper(letter);
    result.insert(result.begin(), letter);
    l /= lmax_plus_1;
  } while (l != 0);

  return result;
}
}  // namespace

CGF::CGF() : pure_sh_(false), unit_(false) {
  for (int i = 0; i < 3; i++) qn_[i] = 0;
}

CGF::CGF(unsigned int qn[3], bool puresh) : pure_sh_(puresh), unit_(false) {
  for (int i = 0; i < 3; i++) qn_[i] = qn[i];
}

CGF::CGF(const CGF& source)
    : Contractable<CGF>(source),
      deriv_(source.deriv_),
      pure_sh_(source.pure_sh_),
      unit_(source.unit_) {
  for (int i = 0; i < 3; i++) qn_[i] = source.qn_[i];
}

CGF::CGF(const ConstructablePolymorphically& sptr)
    : Contractable<CGF>(dynamic_cast<const CGF&>(sptr)) {
  const CGF& sptr_cast = dynamic_cast<const CGF&>(sptr);
  for (int i = 0; i < 3; i++) qn_[i] = sptr_cast.qn_[i];
  deriv_ = sptr_cast.deriv_;
  pure_sh_ = sptr_cast.pure_sh_;
  unit_ = sptr_cast.unit_;
}

CGF::~CGF() {}

std::string CGF::label() const {
  // unit *functions* are treated as regular s functions so that (ss|ss)^(m) =
  // (unit s|ss)^(m)
  // if (is_unit()) return "unit_";
  unsigned int am = qn_[0] + qn_[1] + qn_[2];
  std::string deriv_label;
  if (!deriv_.zero()) deriv_label = deriv_.label();
  const std::string am_string = am_to_symbol(am, contracted());
  std::ostringstream oss;
  oss << (pure_sh_ && am > 0 ? "W" : "") << am_string << deriv_label << "_";
  if (am == 0) return oss.str();

  std::string label = oss.str();
  const char xyz_char[][2] = {"x", "y", "z"};
  for (unsigned int xyz = 0; xyz < 3u; xyz++) {
    std::string xyzlab(xyz_char[xyz]);
    for (unsigned int i = 0; i < qn_[xyz]; i++) {
      label += xyzlab;
    }
  }
  return label;
}

unsigned int CGF::qn(unsigned int i) const {
  assert(i < 3);
  return qn_[i];
}

bool CGF::operator==(const CGF& a) const {
  return (qn_[0] == a.qn_[0] && qn_[1] == a.qn_[1] && qn_[2] == a.qn_[2] &&
          contracted() == a.contracted() && deriv_ == a.deriv_ &&
          pure_sh_ == a.pure_sh_ && unit_ == a.unit_);
}

CGF& CGF::operator=(const CGF& source) {
  for (int i = 0; i < 3; i++) qn_[i] = source.qn_[i];
  deriv_ = source.deriv_;
  pure_sh_ = source.pure_sh_;
  unit_ = source.unit_;
  Contractable<CGF>::operator=(source);
  if (!source.valid()) invalidate();
  return *this;
}

void CGF::dec(unsigned int i, unsigned int c) {
  if (is_unit()) {
    invalidate();
    return;
  }
  if (i < 3 && valid()) {
    if (qn_[i] < c) {
      invalidate();
      return;
    }
    qn_[i] -= c;
  }
}

void CGF::inc(unsigned int i, unsigned int c) {
  assert(!is_unit());
  if (i < 3 && valid()) qn_[i] += c;
}

unsigned int CGF::norm() const { return qn_[0] + qn_[1] + qn_[2]; }

void CGF::print(std::ostream& os) const { os << "CGF: " << label() << endl; }

CGF libint2::operator+(const CGF& A, const CGF& B) {
  assert(!A.is_unit() && !B.is_unit());
  CGF Sum(A);
  for (unsigned int xyz = 0; xyz < 3; ++xyz) Sum.inc(xyz, B.qn(xyz));
  Sum.deriv_ += B.deriv_;
  return Sum;
}

CGF libint2::operator-(const CGF& A, const CGF& B) {
  // assert(A.is_unit() == false && B.is_unit() == false);
  CGF Diff(A);
  for (unsigned int xyz = 0; xyz < 3; ++xyz) Diff.dec(xyz, B.qn(xyz));
  Diff.deriv_ -= B.deriv_;

  return Diff;
}

CGF CGF::unit() {
  CGF result;
  result.unit_ = true;
  result.uncontract();
  return result;
}

///////////////////////////////////////

// By default make it an s-shell
CGShell::CGShell() : pure_sh_(false), unit_(false) {
  for (int i = 0; i < 1; i++) qn_[i] = 0;
}

CGShell::CGShell(unsigned int qn, bool puresh)
    : pure_sh_(puresh), unit_(false) {
  qn_[0] = qn;
}

CGShell::CGShell(const CGShell& source)
    : Contractable<CGShell>(source),
      deriv_(source.deriv_),
      pure_sh_(source.pure_sh_),
      unit_(source.unit_) {
  qn_[0] = source.qn_[0];
}

CGShell::~CGShell() {}

std::string CGShell::label() const {
  if (is_unit()) return "unit";
  std::string result = std::string(pure_sh_ && qn_[0] > 0 ? "W" : "") +
                       am_to_symbol(qn_[0], contracted());
  if (!deriv_.zero()) result += deriv_.label();
  return result;
}

CGShell& CGShell::operator=(const CGShell& source) {
  qn_[0] = source.qn_[0];
  deriv_ = source.deriv_;
  pure_sh_ = source.pure_sh_;
  unit_ = source.unit_;
  Contractable<CGShell>::operator=(source);
  if (!source.valid()) invalidate();
  return *this;
}

bool CGShell::operator==(const CGShell& a) const {
  return (qn_[0] == a.qn_[0] && contracted() == a.contracted() &&
          deriv_ == a.deriv_ && pure_sh_ == a.pure_sh_ && unit_ == a.unit_);
}

void CGShell::dec(unsigned int i, unsigned int c) {
  if (is_unit()) {
    invalidate();
    return;
  }
  if (i == 0 && valid()) {
    if (qn_[0] < c) {
      invalidate();
      return;
    }
    qn_[0] -= c;
  }
}

void CGShell::inc(unsigned int i, unsigned int c) {
  assert(!is_unit());
  if (i == 0 && valid()) qn_[0] += c;
}

unsigned int CGShell::norm() const { return qn_[0]; }

void CGShell::print(std::ostream& os) const {
  os << "CGShell: " << label() << endl;
}

CGShell libint2::operator+(const CGShell& A, const CGShell& B) {
  assert(!A.is_unit() && !B.is_unit());
  CGShell Sum(A);
  Sum.inc(0, B.qn(0));
  Sum.deriv_ += B.deriv_;
  return Sum;
}

CGShell libint2::operator-(const CGShell& A, const CGShell& B) {
  // assert(A.is_unit() == false && B.is_unit() == false);
  CGShell Diff(A);
  Diff.dec(0, B.qn(0));
  Diff.deriv_ -= B.deriv_;
  return Diff;
}

CGShell CGShell::unit() {
  CGShell result;
  result.unit_ = true;
  result.uncontract();
  return result;
}
