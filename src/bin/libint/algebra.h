
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
        OT_(OT), left_(left), right_(right)
        {
        }
      virtual ~AlgebraicOperator() {}
      
      /// Returns the left argument
      SafePtr<T> left() const { return left_; }
      /// Returns the right argument
      SafePtr<T> right() const { return right_; }

      /// Implements DGVertex::size()
      const unsigned int size() const { return 0; }
      /// Implements DGVertex::equiv()
      bool equiv(const SafePtr<DGVertex>& a) const
      {
        SafePtr<AlgebraicOperator> a_cast = dynamic_pointer_cast<AlgebraicOperator,DGVertex>(a);
        if (a_cast) {
          bool result = (OT_ == a_cast->OT_ && left_->equiv(a_cast->left()) && right_->equiv(a_cast->right()));
          return result;
        }
	else
	  return false;
      }
      /// Implements DGVertex::label()
      std::string label() const
      {
        return algebra::OperatorSymbol[OT_];
      }
      /// Implements DGVertex::id()
      std::string id() const
      {
        return label();
      }
      /// Implements DGVertex::description()
      std::string description() const
      {
        ostringstream os;
        os << "AlgebraicOperator: " << left_->description() << " "
           << algebra::OperatorSymbol[OT_] << " "
           << right_->description();
        return os.str();
      }
      /// Implements DGVertex::precomputed
      bool precomputed() const
      {
        return false;
      }

      /// Overloads DGVertex::del_exit_arc()
      void del_exit_arc(const SafePtr<DGArc>& arc)
      {
        throw CannotAddArc("AlgebraicOperator::del_exit_arc() -- cannot safely use del_exit_arc on operator vertices. Use replace_exit_arc instead.");
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
      
    };
    
};

#endif

