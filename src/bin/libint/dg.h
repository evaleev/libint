


#ifndef _libint2_src_bin_libint_dg_h_
#define _libint2_src_bin_libint_dg_h_

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <rr.h>

using namespace std;


namespace libint2 {

  /** DirectedGraph is an implementation of a directed graph
      composed of vertices represented by DGVertex objects. The objects
      are allocated on free store and the graph is implemented as
      vector<DGVertex*>.
  */
  
  class DirectedGraph {

    vector<DGVector*> stack_;

    static const unsigned int default_size_ = 10000;
    unsigned int first_free_;

  public:
    DirectedGraph();
    ~DirectedGraph();

    /** I must derive from DGVertex. RR has a constructor which takes
        const I& as the only argument.

        NOTE TO SELF : need to implement these restrictions using
        standard Bjarne Stroustrup's approach.

    */
    template <class I, class RR> void recurse(const I*);

  }


}


#endif
