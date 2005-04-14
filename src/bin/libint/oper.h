
#include <string>

#ifndef _libint2_src_bin_libint_oper_h_
#define _libint2_src_bin_libint_oper_h_

namespace libint2 {
  
  /** OperatorProperties describes various properties of an operator or operator set
      np -- number of particles
      multi -- true if multiplicative
    */
  template <unsigned int NP, bool multi>
    class OperatorProperties {
    public:
      static const unsigned int np = NP;
      static const bool multiplicative = multi;
    };

  /** OperSet is the base class for all (sets of) operators.
     OperSet's must be constructable using
     SafePtr<OperSet> or SafePtr<ConstructablePolymorphically>.
  */
  class OperSet : public ConstructablePolymorphically {
    public:
      virtual ~OperSet() {};

      /// Returns full description of the operator
      virtual const std::string& descr() const =0;
      /// Returns short label for the operator
      virtual const std::string& label() const =0;

      /** Returns 1, 0, or -1, if each operator in the set is symmetric, nonsymmetric,
        or antisymmetric with respect to permutation of particles i and j */
      virtual int psymm(int i, int j) const =0;

      /// Number of operators in the set
      virtual const unsigned int num_oper() const =0;
    };
  
  /**
  */
  template <class Props>
    class Oper : public OperSet {
    public:
      typedef Props Properties;
      virtual ~Oper();

      /// Returns full description of the operator
      const std::string& descr() const;
      /// Returns short label for the operator
      const std::string& label() const;
      
      bool operator==(const Oper&) const;

    protected:
      /// The only declared constructor is only useable by derived classes
      Oper(const std::string& descr, const std::string& label);

    private:
      /// Described name
      const std::string descr_;
      /// short (<20 chars) ID label
      const std::string label_;
  };

  template <class Props>
    Oper<Props>::Oper(const std::string& descr, const std::string& label) :
    OperSet(), descr_(descr), label_(label)
    {
    }
  
  template <class Props>
    Oper<Props>::~Oper()
    {
    }
  
  template <class Props>
    const std::string&
    Oper<Props>::descr() const
    {
      return descr_;
    }
  
  template <class Props>
    const std::string&
    Oper<Props>::label() const
    {
      return label_;
    }
    
  template <class Props>
  bool
    Oper<Props>::operator==(const Oper& a) const
    {
      return true;
    }


//////////////////////////////
  
  typedef OperatorProperties<2,true> TwoPRep_Props;
  class TwoERep : public Oper<TwoPRep_Props> {
  public:
    typedef Oper<TwoPRep_Props> parent_type;
    /// TwoERep is not a set
    typedef TwoERep iter_type;
    const unsigned int num_oper() const { return 1; };
  
    TwoERep();
    TwoERep(const SafePtr<TwoERep>&);
    TwoERep(const SafePtr<OperSet>&);
    TwoERep(const SafePtr<ConstructablePolymorphically>&);
    ~TwoERep();

    /** Returns 1, 0, or -1, if the operator is symmetric, nonsymmetric,
        or antisymmetric with respect to permutation of particles i and j */
    int psymm(int i, int j) const;

  private:
    // symmetry W.R.T. permutation of each pair of particles
    // 1 -- symmetric, -1 -- antisymmetric, 0 -- nonsymmetric
    // stored as a lower triangle (diagonal not included)
    static const char psymm_[Properties::np*(Properties::np-1)/2];

  };


};

#endif

