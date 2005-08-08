
#include <stdexcept>
#include <cmath>
#include <rr.h>
#include <oper.h>

using namespace std;
using namespace libint2;

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

////////////

std::string
R12kG12::label(int K)
{
  ostringstream oss;
  oss << "R12^" << K << " * G12";
  return oss.str();
}

std::string
R12kG12::symbol(int K)
{
  ostringstream oss;
  oss << "R12_" << (K<0 ? "minus_" : "") << std::abs(K) << "_G12";
  return oss.str();
}

