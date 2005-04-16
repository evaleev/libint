
#include <stdexcept>
#include <smart_ptr.h>

#ifndef _libint2_src_bin_libint_exception_h_
#define _libint2_src_bin_libint_exception_h_


namespace libint2 {

  class DGVertex;

  class InvalidDecrement : public std::logic_error {
    
  public:
    InvalidDecrement(const std::string& a) :
      logic_error(a) {};
    
  };

  class CannotAddArc : public std::logic_error {
    
    public:
    CannotAddArc(const std::string& a) :
      logic_error(a) {};
    
  };

  /** This exception class is used to pass the pointer to the vertex on the graph
   */
  class VertexAlreadyOnStack : public std::logic_error {
    
  public:
    VertexAlreadyOnStack(const SafePtr<DGVertex>& vertex) :
      logic_error("DirectedGraph -- vertex already on stack"), vertex_(vertex) {}
    ~VertexAlreadyOnStack() throw() {}

    SafePtr<DGVertex> vertex() const { return vertex_; }

  private:
    // Vertex on the stack
    SafePtr<DGVertex> vertex_;
    
  };

  /** This exception class is used to notify that a graph operation cannot be performed
   */
  class CannotPerformOperation : public std::logic_error {
  public:
    CannotPerformOperation(const std::string& msg) :
      logic_error(msg) {}
    ~CannotPerformOperation() throw() {}
  };

  /// This exception used to indicate that some property is not set
  template <class T>
    class NotSet : public std::logic_error {

    public:
    NotSet(const std::string& a) :
      logic_error(a) {};
  };


};

#endif

