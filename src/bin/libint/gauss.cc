
#include <bfset.h>
#include <stdexcept>
#include <exception.h>
#include <default_params.h>

using namespace std;
using namespace libint2;

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

CGF::CGF(const CGF& source)
{
  for(int i=0; i<3; i++)
    qn_[i] = source.qn_[i];
}

CGF::CGF(const SafePtr<CGF>& source)
{
  for(int i=0; i<3; i++)
    qn_[i] = source->qn_[i];
}

CGF::CGF(const SafePtr<parent_type>& sptr)
{
  const SafePtr<CGF> sptr_cast = dynamic_pointer_cast<CGF,parent_type>(sptr);
  if (sptr_cast == 0)
    throw std::runtime_error("CGF::CGF(const SafePtr<parent_type>& sptr) -- type of sptr is incompatible with CGF");

  for(int i=0; i<3; i++)
    qn_[i] = sptr_cast->qn_[i];
}

CGF::CGF(const SafePtr<ConstructablePolymorphically>& sptr)
{
  const SafePtr<CGF> sptr_cast = dynamic_pointer_cast<CGF,ConstructablePolymorphically>(sptr);
  if (sptr_cast == 0)
    throw std::runtime_error("CGF::CGF(const SafePtr<ConstructablePolymorphically>& sptr) -- type of sptr is incompatible with CGF");

  for(int i=0; i<3; i++)
    qn_[i] = sptr_cast->qn_[i];
}

CGF::~CGF()
{
}

const std::string
CGF::label() const
{
  unsigned int am = qn_[0] + qn_[1] + qn_[2];
  const char am_char = StaticDefinitions::am_letters[am];
  char tmp[80]; sprintf(tmp,"%c_",am_char);
  if (am == 0) {
    tmp[1] = '\0';
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
           qn_[2] == a.qn_[2] );
}

void
CGF::dec(unsigned int i)
{
  if (i<3) {
    if (qn_[i] == 0)
      throw InvalidDecrement("CGF::dec() -- quantum number already 0");
    --qn_[i];
  }
}

void
CGF::inc(unsigned int i) throw()
{
  if (i<3)
    ++qn_[i];
}

bool
CGF::zero() const
{
  return qn_[0] == 0 && qn_[1] == 0 && qn_[2] == 0;
}

void
CGF::print(std::ostream& os) const
{
  os << "CGF: " << label() << endl;
}

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

CGShell::CGShell(unsigned int qn[1])
{
  for(int i=0; i<1; i++)
    qn_[i] = qn[i];
}

CGShell::CGShell(const CGShell& source)
{
  for(int i=0; i<1; i++)
    qn_[i] = source.qn_[i];
}

CGShell::CGShell(const SafePtr<CGShell>& source)
{
  for(int i=0; i<1; i++)
    qn_[i] = source->qn_[i];
}

CGShell::CGShell(const SafePtr<parent_type>& sptr)
{
  const SafePtr<CGShell> sptr_cast = dynamic_pointer_cast<CGShell,parent_type>(sptr);
  if (sptr_cast == 0)
    throw std::runtime_error("CGShell::CGShell(const SafePtr<parent_type>& sptr) -- type of sptr is incompatible with CGShell");

  for(int i=0; i<1; i++)
    qn_[i] = sptr_cast->qn_[i];
}

CGShell::CGShell(const SafePtr<ConstructablePolymorphically>& sptr)
{
  const SafePtr<CGShell> sptr_cast = dynamic_pointer_cast<CGShell,ConstructablePolymorphically>(sptr);
  if (sptr_cast == 0)
    throw std::runtime_error("CGShell::CGShell(const SafePtr<ConstructablePolymorphically>& sptr) -- type of sptr is incompatible with CGShell");

  for(int i=0; i<1; i++)
    qn_[i] = sptr_cast->qn_[i];
}

CGShell::~CGShell()
{
}

const std::string
CGShell::label() const
{
  return std::string(1,StaticDefinitions::am_letters[qn_[0]]);
}

CGShell&
CGShell::operator=(const CGShell& source)
{
  for(int i=0; i<1; i++)
    qn_[i] = source.qn_[i];
  return *this;
}

bool
CGShell::operator==(const CGShell& a) const
{
  return ( qn_[0] == a.qn_[0] );
}

void
CGShell::dec(unsigned int i)
{
  if (i == 0) {
    if (qn_[0] == 0)
      throw InvalidDecrement("CGShell::dec() -- angular momentum already 0");
    --qn_[0];
  }
  else {
    throw InvalidDecrement("CGShell::dec(i) -- i must be zero");
  }
}

void
CGShell::inc(unsigned int i) throw()
{
  if (i == 0)
    ++qn_[0];
}

bool
CGShell::zero() const
{
  return qn_[0] == 0;
}

void
CGShell::print(std::ostream& os) const
{
  os << "CGShell: am = " << qn_[0] << endl;
}

