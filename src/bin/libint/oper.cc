
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

////////////

R1dotR1_G12::R1dotR1_G12() :
  parent_type("r_1.r_1 x G12","R1dotR1_G12")
{
}

R1dotR1_G12::R1dotR1_G12(const SafePtr<R1dotR1_G12>& source) :
  parent_type("r_1.r_1 x G12","R1dotR1_G12")
{
}

R1dotR1_G12::R1dotR1_G12(const SafePtr<OperSet>& oset) :
  parent_type("r_1.r_1 x G12","R1dotR1_G12")
{
  const SafePtr<R1dotR1_G12> oset_cast = dynamic_pointer_cast<R1dotR1_G12,OperSet>(oset);
  if (oset_cast == 0)
    throw std::runtime_error("R1dotR1_G12::R1dotR1_G12(const SafePtr<OperSet>& oset) -- oset is a pointer to an incompatible type");
}

R1dotR1_G12::R1dotR1_G12(const SafePtr<ConstructablePolymorphically>& oset) :
  parent_type("r_1.r_1 x G12","R1dotR1_G12")
{
  const SafePtr<R1dotR1_G12> oset_cast = dynamic_pointer_cast<R1dotR1_G12,ConstructablePolymorphically>(oset);
  if (oset_cast == 0)
    throw std::runtime_error("R1dotR1_G12::R1dotR1_G12(const SafePtr<ConstructablePolymorphically>& oset) -- oset is a pointer to an incompatible type");
}

R1dotR1_G12::~R1dotR1_G12()
{
}

////////////

R1dotR2_G12::R1dotR2_G12() :
  parent_type("r_1.r_2 x G12","R1dotR2_G12")
{
}

R1dotR2_G12::R1dotR2_G12(const SafePtr<R1dotR2_G12>& source) :
  parent_type("r_1.r_2 x G12","R1dotR2_G12")
{
}

R1dotR2_G12::R1dotR2_G12(const SafePtr<OperSet>& oset) :
  parent_type("r_1.r_2 x G12","R1dotR2_G12")
{
  const SafePtr<R1dotR2_G12> oset_cast = dynamic_pointer_cast<R1dotR2_G12,OperSet>(oset);
  if (oset_cast == 0)
    throw std::runtime_error("R1dotR2_G12::R1dotR2_G12(const SafePtr<OperSet>& oset) -- oset is a pointer to an incompatible type");
}

R1dotR2_G12::R1dotR2_G12(const SafePtr<ConstructablePolymorphically>& oset) :
  parent_type("r_1.r_2 x G12","R1dotR2_G12")
{
  const SafePtr<R1dotR2_G12> oset_cast = dynamic_pointer_cast<R1dotR2_G12,ConstructablePolymorphically>(oset);
  if (oset_cast == 0)
    throw std::runtime_error("R1dotR2_G12::R1dotR2_G12(const SafePtr<ConstructablePolymorphically>& oset) -- oset is a pointer to an incompatible type");
}

R1dotR2_G12::~R1dotR2_G12()
{
}

////////////
