
#include <stdexcept>
#include <rr.h>

using namespace std;
using namespace libint2;

const char TwoERep::psymm_[] = { 1 };

TwoERep::TwoERep() :
  Oper<TwoPRep_Props>("Two-electron repulsion operator","TwoERep")
{
}

TwoERep::~TwoERep()
{
}

const int
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

