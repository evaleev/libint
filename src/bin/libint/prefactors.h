
#include <smart_ptr.h>
#include <entity.h>

#ifndef _libint2_src_bin_libint_prefactors_h_
#define _libint2_src_bin_libint_prefactors_h_

namespace libint2 {
  
  /** 
     Prefactors is a collection of common quantities which appear
     as prefactors in recurrence relations for Gaussian integrals.
  */
  
  class Prefactors {
    
    public:
    Prefactors();
    ~Prefactors();
    
    typedef RTimeEntity<double> rdouble;

    /// P-A vector
    SafePtr<rdouble> PA;
    /// Q-C vector
    SafePtr<rdouble> QC;
    /// W-P vector
    SafePtr<rdouble> WP;
    /// W-Q vector
    SafePtr<rdouble> WQ;

  };
  
  extern Prefactors prefactors;
  
};

#endif

