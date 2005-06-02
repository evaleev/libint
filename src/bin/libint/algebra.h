
#include <smart_ptr.h>
#include <rr.h>
#include <exception.h>

#ifndef _libint2_src_bin_libint_algebra_h_
#define _libint2_src_bin_libint_algebra_h_

namespace libint2 {

  namespace algebra {
    struct OperatorTypes {
      typedef enum {Plus, Minus, Times, Divide} OperatorType;
    };
    static const char OperatorSymbol[][2] = { "+", "-", "*", "/" };
  };


  /**
    AlgebraicOperator is an algebraic operator that
    acts on objects of type T. An algebraic operator has 2 arguments,
    to the left and to the right ( L oper R ).
  */

  template <class T>
    class AlgebraicOperator :
    public DGVertex
    {
      
    public:
      typedef algebra::OperatorTypes OperatorTypes;
      typedef algebra::OperatorTypes::OperatorType OperatorType;

      AlgebraicOperator(OperatorType OT,
                        const SafePtr<T>& left,
                        const SafePtr<T>& right) :
        DGVertex(ClassInfo<AlgebraicOperator>::Instance().id()), OT_(OT), left_(left), right_(right),
        descr_(), label_(algebra::OperatorSymbol[OT_])
        {
        }
      virtual ~AlgebraicOperator() {}
      
      /// Returns the left argument
      SafePtr<T> left() const { return left_; }
      /// Returns the right argument
      SafePtr<T> right() const { return right_; }

      /// Implements DGVertex::size()
      const unsigned int size() const { return 1; }
      /// Implements DGVertex::equiv()
      bool equiv(const SafePtr<DGVertex>& a) const
      {
        if (typeid_ == a->typeid_) {
          // Safe to cast statically because type of a has just been checked
          SafePtr<AlgebraicOperator> a_cast = static_pointer_cast<AlgebraicOperator,DGVertex>(a);
          if (OT_ == a_cast->OT_ && left_->equiv(a_cast->left()) && right_->equiv(a_cast->right()))
            return true;
          else
            return false;
        }
	else
	  return false;
      }
      /// Implements DGVertex::label()
      const std::string& label() const
      {
        return label_;
      }
      /// Implements DGVertex::id()
      const std::string& id() const
      {
        return label();
      }
      /// Implements DGVertex::description()
      const std::string& description() const
      {
        if (descr_.empty()) {
          ostringstream os;
          os << "AlgebraicOperator: " << left_->description() << " "
             << algebra::OperatorSymbol[OT_] << " "
             << right_->description();
          descr_ = os.str();
        }
        return descr_;
      }

      /// Overloads DGVertex::del_exit_arcs()
      void del_exit_arcs()
      {
        throw CannotAddArc("AlgebraicOperator::del_exit_arcs() -- cannot safely use del_exit_arcs on operator vertices.");
      }
      
    private:
      OperatorType OT_;
      SafePtr<T> left_;
      SafePtr<T> right_;

      /// Implements DGVertex::this_precomputed()
      bool this_precomputed() const
      {
        return false;
      }

      std::string label_;
      mutable std::string descr_;
    };
    
};

#endif

