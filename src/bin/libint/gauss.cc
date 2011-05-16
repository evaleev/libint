
#include <ctype.h>
#include <bfset.h>
#include <stdexcept>
#include <exception.h>
#include <default_params.h>

using namespace std;
using namespace libint2;

unsigned CGF::key_l_offset[] = { 0, 1, 4, 10, 20, 35, 56, 84, 120, 165, 220, 286, 364, 455, 560, 680, 816, 969, 1140, 1330, 1540};

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

CGF::CGF(const CGF& source) : Contractable<CGF>(source)
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
}

CGF::~CGF()
{
}

const std::string
CGF::label() const
{
  unsigned int am = qn_[0] + qn_[1] + qn_[2];
  const char am_char = StaticDefinitions::am_letters[am];
  char tmp[80]; sprintf(tmp,"%c_", contracted() ? toupper(am_char) : am_char);
  // to differentiate s-type CGF from s-type CGShell, use "s_"
  if (am == 0) {
    tmp[2] = '\0';
    return tmp;
  }
  std::string label(tmp);
  const char xyz_char[][2] = {"x","y","z"};
  for(int xyz=0; xyz<3; xyz++) {
    std::string xyzlab(xyz_char[xyz]);
    for(int i=0; i<qn_[xyz]; i++) {
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
           contracted() == a.contracted());
}

CGF&
CGF::operator=(const CGF& source)
{
  for(int i=0; i<3; i++)
    qn_[i] = source.qn_[i];
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

CGShell::CGShell(const CGShell& source) : Contractable<CGShell>(source)
{
  for(int i=0; i<1; i++)
    qn_[i] = source.qn_[i];
}

CGShell::~CGShell()
{
}

const std::string
CGShell::label() const
{
  const char amchar = StaticDefinitions::am_letters[qn_[0]];
  return std::string(1,contracted() ? toupper(amchar) : amchar);
}

CGShell&
CGShell::operator=(const CGShell& source)
{
  for(int i=0; i<1; i++)
    qn_[i] = source.qn_[i];
  Contractable<CGShell>::operator=(source);
  if (!source.valid()) invalidate();
  return *this;
}

bool
CGShell::operator==(const CGShell& a) const
{
  return ( qn_[0] == a.qn_[0] &&
           contracted() == a.contracted() );
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
  return Sum;
}

CGShell
libint2::operator-(const CGShell& A, const CGShell& B) {
  CGShell Diff(A);
  Diff.dec(0,B.qn(0));
  return Diff;
}
