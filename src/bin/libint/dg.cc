
#include <rr.h>

using namespace std;
using namespace libint2;

DGArc::DGArc(const DGVertex* orig, const DGVertex* dest)
{
  orig_ = orig;
  dest_ = dest;
}

DGArc::~DGArc()
{
}

DGVertex::DGVertex() :
  parents_(0), children_(0)
{
}

DGVertex::DGVertex(const vector<DGArc*>& parents, const vector<DGArc*>& children) :
  parents_(parents), children_(children)
{
}

DGVertex::~DGVertex()
{
}

