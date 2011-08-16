
#ifndef _libint2_src_bin_libint_key_h_
#define _libint2_src_bin_libint_key_h_

#include <libint2_intrinsic_types.h>
#include <gmpxx.h>

namespace libint2 {

  /// Type/Instance combination serves as a key to quickly compare 2 polymorphic Singletons
  template <typename T, typename I>
  class TypeAndInstance {
  public:
    typedef T Type;
    typedef I Instance;
    TypeAndInstance() : t_(invalid_type_), i_(invalid_instance_) {}
    TypeAndInstance(const Type& t, const Instance& i) : t_(t), i_(i) {}
    TypeAndInstance(const TypeAndInstance& i) : t_(i.t_), i_(i.i_) {}
    const TypeAndInstance& operator=(const TypeAndInstance& i) { t_ = i.t_; i_ = i.i_; return *this; }

    const Type& type() const { return t_; }
    const Instance& instance() const { return i_; }
    
  private:
    Type t_;
    Instance i_;
    
    static Type invalid_type_;
    static Instance invalid_instance_;
  };

  template <typename T, typename I> typename TypeAndInstance<T,I>::Type TypeAndInstance<T,I>::invalid_type_(-1);
  template <typename T, typename I> typename TypeAndInstance<T,I>::Instance TypeAndInstance<T,I>::invalid_instance_(-1);

  template <typename T, typename I>
  bool operator==(const TypeAndInstance<T,I>& a,
		  const TypeAndInstance<T,I>& b) {
    return a.type() == b.type() && a.instance() == b.instance();
  }
  
  template <typename T, typename I>
  bool operator<(const TypeAndInstance<T,I>& a,
		 const TypeAndInstance<T,I>& b) {
    bool result = 
      (a.type() < b.type()) || 
      ( (a.type() == b.type()) &&
	(a.instance() < b.instance())
      );
    return result;
  }


  ///
  /// Collection of types used for constructing keys in libint2
  ///
  struct KeyTypes {
    /// distinct classes have unique ClassID's
    typedef unsigned int ClassID;
    /// some classes need to have distinct instances to have unique InstanceID's, e.g. generalized Singletons
    typedef mpz_class InstanceID;
  };
  /// this composite hashing key works for DGVertex
  typedef TypeAndInstance<KeyTypes::ClassID,KeyTypes::InstanceID> DGVertexKey;

};

#endif

// Local Variables:
// mode: c++
// End:
