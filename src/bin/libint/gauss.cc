
#include <rr.h>

using namespace std;
using namespace libint2;

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

CGF::CGF(const BFSet* sptr)
{
  const CGF* sptr_cast = dynamic_cast<const CGF*>(sptr);
  if (sptr_cast == 0)
    throw std::runtime_error("CGF::CGF(const BFSet* sptr) -- type of sptr is incompatible with CGF");

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
  char* const c_label = new char[80];
  sprintf(c_label,"%c(x^%d y^%d z^%d)",am_char,qn_[0],qn_[1],qn_[2]);
  return c_label;
}

bool
CGF::operator==(const CGF& a) const
{
  return ( qn_[0] == a.qn_[0] &&
           qn_[1] == a.qn_[1] &&
           qn_[2] == a.qn_[2] );
}


// By default make it an s-shell
CGShell::CGShell()
{
  for(int i=0; i<1; i++)
    qn_[i] = 0;
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

CGShell::CGShell(const BFSet* sptr)
{
  const CGShell* sptr_cast = dynamic_cast<const CGShell*>(sptr);
  if (sptr_cast == 0)
    throw std::runtime_error("CGShell::CGShell(const BFSet* sptr) -- type of sptr is incompatible with CGShell");

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
}

bool
CGShell::operator==(const CGShell& a) const
{
  return ( qn_[0] == a.qn_[0] );
}

void
CGShell::dec()
{
  // NOTE TO SELF : need to throw specialized exception
  if (qn_[0] == 0)
    throw InvalidDecrement("CGShell::dec() -- angular momentum already 0");
  --qn_[0];
}

void
CGShell::inc() throw()
{
  ++qn_[0];
}

void
CGShell::print(std::ostream& os) const
{
  os << "CGShell: am = " << qn_[0] << endl;
}
