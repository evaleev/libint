
#include <iostream>
#include <sstream>
#include <ctype.h>
#include <bfset.h>
#include <stdexcept>
#include <string>
#include <exception.h>
#include <default_params.h>
#include <string.h>

using namespace std;
using namespace libint2;

unsigned CGF::key_l_offset[] = { 0, 1, 4, 10, 20, 35, 56, 84, 120, 165, 220, 286, 364, 455, 560, 680, 816, 969, 1140, 1330, 1540,
                                 1771, 2024, 2300, 2600, 2925, 3276, 3654, 4060, 4495};
unsigned OriginDerivative::key_l_offset[] = { 0, 1, 4, 10, 20};

namespace {
  std::string am_to_symbol(unsigned int l, bool contracted) {
    std::string result;
    do {
      const unsigned int digit = l % 10u;
      char letter = StaticDefinitions::am_letters[digit];
      if (contracted)
        letter = toupper(letter);
      result.insert(result.begin(), letter);
      l /= 10;
    } while (l != 0);

    return result;
  }
}

OriginDerivative
libint2::operator-(const OriginDerivative& A, const OriginDerivative& B) {
  OriginDerivative Diff(A);
  for(unsigned int xyz=0; xyz<3; ++xyz)
    Diff.dec(xyz,B.d(xyz));
  return Diff;
}

bool
libint2::operator==(const OriginDerivative& A, const OriginDerivative& B) {
  for(unsigned int xyz=0; xyz<3; ++xyz)
    if (A.d(xyz) != B.d(xyz))
      return false;
  return true;
}

CGF::CGF()
{
  for(int i=0; i<3; i++)
    qn_[i] = 0;
}

CGF::CGF(unsigned int qn[3])
{
  for(int i=0; i<3; i++)
    qn_[i] = qn[i];
}

CGF::CGF(const CGF& source) : Contractable<CGF>(source),
    deriv_(source.deriv_)
{
  for(int i=0; i<3; i++)
    qn_[i] = source.qn_[i];
}

CGF::CGF(const ConstructablePolymorphically& sptr) :
    Contractable<CGF>(dynamic_cast<const CGF&>(sptr))
{
  const CGF& sptr_cast = dynamic_cast<const CGF&>(sptr);
  for(int i=0; i<3; i++)
    qn_[i] = sptr_cast.qn_[i];
  deriv_ = sptr_cast.deriv_;
}

CGF::~CGF()
{
}

const std::string
CGF::label() const
{
  unsigned int am = qn_[0] + qn_[1] + qn_[2];
  std::string deriv_label;
  if (deriv_.zero() == false) deriv_label = deriv_.label();
  const std::string am_string = am_to_symbol(am, contracted());
  std::ostringstream oss;
  oss << am_string << deriv_label << "_";
  if (am == 0) return oss.str();

  std::string label = oss.str();
  const char xyz_char[][2] = {"x","y","z"};
  for(unsigned int xyz=0; xyz<3u; xyz++) {
    std::string xyzlab(xyz_char[xyz]);
    for(unsigned int i=0; i<qn_[xyz]; i++) {
      label += xyzlab;
    }
  }
  return label;
}

unsigned int
CGF::qn(unsigned int i) const
{
  assert(i < 3);
  return qn_[i];
}

bool
CGF::operator==(const CGF& a) const
{
  return ( qn_[0] == a.qn_[0] &&
           qn_[1] == a.qn_[1] &&
           qn_[2] == a.qn_[2] &&
           contracted() == a.contracted() &&
           deriv_ == a.deriv_);
}

CGF&
CGF::operator=(const CGF& source)
{
  for(int i=0; i<3; i++)
    qn_[i] = source.qn_[i];
  deriv_ = source.deriv_;
  Contractable<CGF>::operator=(source);
  if (!source.valid()) invalidate();
  return *this;
}

void
CGF::dec(unsigned int i, unsigned int c)
{
  if (i<3 && valid()) {
    if (qn_[i] < c) {
      invalidate();
      return;
    }
    qn_[i] -= c;
  }
}

void
CGF::inc(unsigned int i, unsigned int c)
{
  if (i<3 && valid())
    qn_[i] += c;
}

unsigned int
CGF::norm() const
{
  return qn_[0] + qn_[1] + qn_[2];
}

void
CGF::print(std::ostream& os) const
{
  os << "CGF: " << label() << endl;
}

CGF
libint2::operator+(const CGF& A, const CGF& B) {
  CGF Sum(A);
  for(unsigned int xyz=0; xyz<3; ++xyz)
    Sum.inc(xyz,B.qn(xyz));
  Sum.deriv_ += B.deriv_;
  return Sum;
}

CGF
libint2::operator-(const CGF& A, const CGF& B) {
  CGF Diff(A);
  for(unsigned int xyz=0; xyz<3; ++xyz)
    Diff.dec(xyz,B.qn(xyz));

  return Diff;
}

///////////////////////////////////////

// By default make it an s-shell
CGShell::CGShell()
{
  for(int i=0; i<1; i++)
    qn_[i] = 0;
}

CGShell::CGShell(unsigned int qn)
{
    qn_[0] = qn;
}

CGShell::CGShell(const CGShell& source) : Contractable<CGShell>(source),
    deriv_(source.deriv_)
{
    qn_[0] = source.qn_[0];
}

CGShell::~CGShell()
{
}

const std::string
CGShell::label() const
{
  std::string result = am_to_symbol(qn_[0], contracted());
  if (!deriv_.zero())
    result += deriv_.label();
  return result;
}

CGShell&
CGShell::operator=(const CGShell& source)
{
  qn_[0] = source.qn_[0];
  deriv_ = source.deriv_;
  Contractable<CGShell>::operator=(source);
  if (!source.valid()) invalidate();
  return *this;
}

bool
CGShell::operator==(const CGShell& a) const
{
  return ( qn_[0] == a.qn_[0] &&
           contracted() == a.contracted() &&
           deriv_ == a.deriv_);
}

void
CGShell::dec(unsigned int i, unsigned int c)
{
  if (i == 0 && valid()) {
    if (qn_[0] < c) {
      invalidate();
      return;
    }
    qn_[0] -= c;
  }
}

void
CGShell::inc(unsigned int i, unsigned int c)
{
  if (i == 0 && valid())
    qn_[0] += c;
}

unsigned int
CGShell::norm() const
{
  return qn_[0];
}

void
CGShell::print(std::ostream& os) const
{
  os << "CGShell: " << label() << endl;
}

CGShell
libint2::operator+(const CGShell& A, const CGShell& B) {
  CGShell Sum(A);
  Sum.inc(0,B.qn(0));
  Sum.deriv_ += B.deriv_;
  return Sum;
}

CGShell
libint2::operator-(const CGShell& A, const CGShell& B) {
  CGShell Diff(A);
  Diff.dec(0,B.qn(0));
  Diff.deriv_ -= B.deriv_;
  return Diff;
}
