
#include <smart_ptr.h>
#include <rr.h>

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
          cout << "AlgebraicOperator::equiv:" << endl << "  ";
          print(cout); cout << endl << "  ";
          a_cast->print(cout); cout << endl;
          cout << "  " << (result ? "true" : "false") << endl;
          return result;
        }
	else
	  return false;
      }
      /// Implements DGVertex::print()
      void print(std::ostream& os) const
      {
	/*os << "AlgebraicOperator: ";
        left_->print(os);
        os << " " << algebra::OperatorSymbol[OT_] << " ";
        right_->print(os);*/
        os <<  algebra::OperatorSymbol[OT_];
      }
      /// Implements DGVertex::precomputed
      bool precomputed() const
      {
        return false;
      }
      
    private:
      OperatorType OT_;
      SafePtr<T> left_;
      SafePtr<T> right_;
      
    };
    
};

#endif

