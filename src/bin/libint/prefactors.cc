
#include <smart_ptr.h>
#include <prefactors.h>

using namespace libint2;

Prefactors::Prefactors() :
  PA(new rdouble("PA")),
  QC(new rdouble("QC")),
  WP(new rdouble("WP")),
  WQ(new rdouble("WQ"))
{
}

Prefactors::~Prefactors()
{
}

Prefactors prefactors;

