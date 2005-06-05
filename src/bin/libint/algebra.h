
#include <smart_ptr.h>
#include <rr.h>
#include <exception.h>
#include <global_macros.h>
#include <dgvertex.h>

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

      /// Clone A but replace operands with left and right
      AlgebraicOperator(const SafePtr<AlgebraicOperator>& A,
                        const SafePtr<T>& left,
                        const SafePtr<T>& right) :
        DGVertex(static_cast<DGVertex&>(*A)), OT_(A->OT_),
        left_(left), right_(right), descr_(), label_(A->label_)
        {
#if DEBUG
          if (num_exit_arcs() != 2)
            cout << "AlgebraicOperator<DGVertex> copy constructor: number of children != 2" << endl;
          else {
            if (left_ != exit_arc(0)->dest() && left_ != exit_arc(1)->dest())
              cout << "AlgebraicOperator<DGVertex> copy constructor: invalid left operand given" << endl; 
            if (right_ != exit_arc(0)->dest() && right_ != exit_arc(1)->dest())
              cout << "AlgebraicOperator<DGVertex> copy constructor: invalid right operand given" << endl;
          }
#endif
        }
      
      /// Returns the left argument
      const SafePtr<T>& left() const { return left_; }
      /// Returns the right argument
      const SafePtr<T>& right() const { return right_; }

      /// Overloads DGVertex::add_exit_arc(). The default definition is used unless T = DGVertex (see algebra.cc)
      void add_exit_arc(const SafePtr<DGArc>& a);
      /// Implements DGVertex::size()
      const unsigned int size() const { return 1; }
      /// Implements DGVertex::equiv()
      bool equiv(const SafePtr<DGVertex>& a) const
      {
        if (typeid_ == a->typeid_) {
#if ALGEBRAICOPERATOR_USE_KEY_TO_COMPARE
          return description() == a->description();
#else
          // Safe to cast statically because type of a has just been checked
          SafePtr<AlgebraicOperator> a_cast = static_pointer_cast<AlgebraicOperator,DGVertex>(a);

#if 1
          // Find out why sometimes equivalent left_ and a_cast->left_ have non-equivalent pointers
          if (left_->equiv(a_cast->left()) && left_ != a_cast->left_) {
            cout << "Left arguments are equivalent but pointers differ!" << endl;
            cout << left_->description() << endl;
            cout << a_cast->left_->description() << endl;
          }
          // Find out why sometimes equivalent right_ and a_cast->right_ have non-equivalent pointers
          if (right_->equiv(a_cast->right()) && right_ != a_cast->right_) {
            cout << "Left arguments are equivalent but pointers differ!" << endl;
            cout << right_->description() << endl;
            cout << a_cast->right_->description() << endl;
          }
#endif
          if (OT_ == a_cast->OT_) {
  #if ALGEBRAICOPERATOR_USE_SAFEPTR
            if (left_ == a_cast->left_ && right_ == a_cast->right_)
  #else
            if (left_->equiv(a_cast->left()) && right_->equiv(a_cast->right()))
  #endif
            return true;
          else
            return false;
          }
          else
            return false;
#endif
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

  /*
  template <>
    void
    AlgebraicOperator<DGVertex>::add_exit_arc(const SafePtr<DGArc>& a)
    {
      DGVertex::add_exit_arc(a);
      if (left_->equiv(a->dest()))
        left_ = a->dest();
      else if (right_->equiv(a->dest()))
        right_ = a->dest();
      else
        throw std::runtime_error("AlgebraicOperator<DGVertex>::add_exit_arc -- trying to add an arc to a vertex not equivalent to either argument.");
    }
  */
  
};

#endif

