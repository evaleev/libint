
#include <map>
#include <dgvertex.h>
#include <entity.h>
#include <extract.h>

using namespace std;
using namespace libint2;

void
ExtractPrecomputedLabels::operator()(const VertexPtr& v)
{
  if (v->precomputed()) {

    // discard compile-time entities
    {
      typedef CTimeEntity<double> cdouble;
      SafePtr<cdouble> ptr_cast = dynamic_pointer_cast<cdouble,DGVertex>(v);
      if (ptr_cast) {
	return;
      }
    }

    map_[v->label()] = true;
  }

}

const ExtractPrecomputedLabels::Labels&
ExtractPrecomputedLabels::labels()
{
  labels_.clear();
  typedef LabelMap::const_iterator citer;
  citer end = map_.end();
  for(citer l=map_.begin(); l!=end; ++l) {
    labels_.push_back(l->first);
  }
  labels_.sort();
  return labels_;
}
