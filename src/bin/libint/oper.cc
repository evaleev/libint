
#include <stdexcept>
#include <rr.h>
#include <oper.h>

using namespace std;
using namespace libint2;

const char TwoERep::psymm_[] = { 1 };

TwoERep::TwoERep() :
  parent_type("Two-electron repulsion operator","TwoERep")
{
}

TwoERep::TwoERep(const SafePtr<TwoERep>& source) :
  parent_type("Two-electron repulsion operator","TwoERep")
{
}

TwoERep::TwoERep(const SafePtr<OperSet>& oset) :
  parent_type("Two-electron repulsion operator","TwoERep")
{
  const SafePtr<TwoERep> oset_cast = dynamic_pointer_cast<TwoERep,OperSet>(oset);
  if (oset_cast == 0)
    throw std::runtime_error("TwoERep::TwoERep(const SafePtr<OperSet>& oset) -- oset is a pointer to an incompatible type");
}

TwoERep::TwoERep(const SafePtr<ConstructablePolymorphically>& oset) :
  parent_type("Two-electron repulsion operator","TwoERep")
{
  const SafePtr<TwoERep> oset_cast = dynamic_pointer_cast<TwoERep,ConstructablePolymorphically>(oset);
  if (oset_cast == 0)
    throw std::runtime_error("TwoERep::TwoERep(const SafePtr<ConstructablePolymorphically>& oset) -- oset is a pointer to an incompatible type");
}

TwoERep::~TwoERep()
{
}

int
TwoERep::psymm(int i, int j) const
{
  if (i<0 || i>=Properties::np)
    throw std::runtime_error("TwoERep::psymm(i,j) -- index i out of bounds");
  if (j<0 || j>=Properties::np)
    throw std::runtime_error("TwoERep::psymm(i,j) -- index j out of bounds");
  if (i == j)
    return 1;
  int ii = (i > j) ? i : j;
  int jj = (i > j) ? j : i;
  return psymm_[ii*(ii-1)/2 + jj];
}

////////////

