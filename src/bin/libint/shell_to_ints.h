
#ifndef _libint2_src_bin_libint_shelltoints_h_
#define _libint2_src_bin_libint_shelltoints_h_

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <rr.h>
#include <integral.h>

using namespace std;


namespace libint2 {

  /** CGShell_to_CGFs converts a shell set to individual integrals. Although this is
  technically not a recurrence relation, it can be treated as one.
  */
  template <template <class> class I>
  class CGShell_to_CGF : public RecurrenceRelation {

    I<CGShell>* target_;
    I<CGF>** children_;

    unsigned int num_actual_children_;

  public:
    CGShell_to_CGF(I<CGShell>*);
    ~CGShell_to_CGF();

    typedef I<CGShell> TargetType;

    const unsigned int num_children() const { return num_actual_children_; };
    /// target() returns points to the i-th child
    I<CGShell>* target() { return target_; };
    /// child(i) returns points i-th child
    I<CGF>* child(unsigned int i);

    const std::string cpp_function_name() {};
    const std::string cpp_source_name() {};
    const std::string cpp_header_name() {};
    std::ostream& cpp_source(std::ostream&) {};

  };
  
  
  template <template <class> class I>
    CGShell_to_CGF<I>::CGShell_to_CGF(I<CGShell>* Tint) :
    target_(Tint)
    {
      target_ = Tint;
      unsigned int m = Tint->m();

      typedef typename I<CGShell>::BraType IBraType;
      typedef typename I<CGShell>::KetType IKetType;
      IBraType* bra = Tint->bra();
      IKetType* ket = Tint->ket();

      // Compute the number of children
      num_actual_children_ = 1;
      typedef typename I<CGShell>::OperatorType OperType;
      unsigned int num_part = OperType::Properties::np;
      for(int p=0; p<num_part; p++) {
        const unsigned int nfunc_bra = bra->num_members(p);
        const unsigned int nfunc_ket = ket->num_members(p);
        for(unsigned int i=0; i<nfunc_bra; i++)
          num_actual_children_ *= bra->member(p,i)->num_bf();
        for(unsigned int i=0; i<nfunc_ket; i++)
          num_actual_children_ *= ket->member(p,i)->num_bf();
      }
      
      // Zero out children pointers
      children_ = new I<CGF>*[num_actual_children_];
      for(int i=0; i<num_actual_children_; i++)
        children_[i] = 0;

      cout << "CGShell_to_CGF<I>::CGShell_to_CGF -- number of children = "
        << num_actual_children_ << endl;
      
      /// use iterator<I> to iterate over elements
    };

  template <template <class> class I>
    I<CGF>*
    CGShell_to_CGF<I>::child(unsigned int i)
  {
      assert(i>=0 && i<num_actual_children_);
      return children_[i];
  };
  

};

#endif
