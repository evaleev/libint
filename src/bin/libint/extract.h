
/**
   Classes here extract information from DirectedGraph
*/

#ifndef _libint2_src_bin_libint_extract_h_
#define _libint2_src_bin_libint_extract_h_

#include <string>
#include <list>
#include <smart_ptr.h>

namespace libint2 {

  class DGVertex;

  /// This class collects labels of all precomputed non-compile-time constants
  class ExtractPrecomputedLabels {
  public:
    typedef SafePtr<DGVertex> VertexPtr;
    typedef std::list<std::string> Labels;

    ExtractPrecomputedLabels() {}
    ~ExtractPrecomputedLabels() {}

    /// try v
    void operator()(const VertexPtr& v);

    /// return list of sorted labels
    const Labels& labels();

  private:
    Labels labels_;
    // labels are stored as a map
    typedef std::map<std::string,bool> LabelMap;
    LabelMap map_;
  };


};

#endif // header guard
