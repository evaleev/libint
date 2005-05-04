
#include <smart_ptr.h>
#include <entity.h>

#ifndef _libint2_src_bin_libint_dims_h_
#define _libint2_src_bin_libint_dims_h_

using namespace std;

namespace libint2 {
  
  using namespace EntityTypes;
  
  /** ImplicitDimensions describes basis functions or other "degrees of freedom"
      not actively engaged in a recurrence relation. For example, horizontal
      AM transfer to obtain a (dd|ps) integral from (fp|ps) and (dp|ps) does
      not involve the |ps) part. Thus there's may be no reason to generate transfer
      routine specific to this integral, only a routine specific to the (dd| part.
      Such function will require the information about the rank of the |ps) part.
      This information is encoded in ImplicitDimensions.
      
      Another special dimension is the vector length...
  */
  
  class ImplicitDimensions {
    public:    
    /// Explicitly initialize both quantities. Their exact time is not known.
    ImplicitDimensions(const SafePtr<Entity>& high,
                       const SafePtr<Entity>& low,
                       const SafePtr<Entity>& vecdim);
    /// Default assumes runtime (dynamical) quantities
    ImplicitDimensions();
    /// Handy constructor to initialize dimensions as compile-time (static) quatities
    ImplicitDimensions(int high, int low, int vec);
    ~ImplicitDimensions() {}
    
    /// Returns the high dimension
    SafePtr<Entity> high() const { return high_; }
    /// Returns the low dimension
    SafePtr<Entity> low() const { return low_; }
    /// Returns the vector dimension
    SafePtr<Entity> vecdim() const { return vecdim_; }
    /// Returns true if the rank of high dimension is known
    bool high_is_static() const { return high_is_static_; }
    /// Returns true if the rank of low dimension is known
    bool low_is_static() const { return low_is_static_; }
    /// Returns true if the rank of vector dimension is known
    bool vecdim_is_static() const { return vecdim_is_static_; }
    
    /// Default ImplicitDimension object
    static SafePtr<ImplicitDimensions> default_dims();
    
    private:
    // Dimensions can be runtime or compile-time quantities
    SafePtr<Entity> high_;
    SafePtr<Entity> low_;
    SafePtr<Entity> vecdim_;
    
    // checks if the dimensions are CTImeEntities
    void init_();
    bool high_is_static_;
    bool low_is_static_;
    bool vecdim_is_static_;
    
    /// Default dimension
    static SafePtr<ImplicitDimensions> default_dims_;
    
  };
  
};

#endif

