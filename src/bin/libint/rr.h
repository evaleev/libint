

#ifndef _libint2_src_bin_libint_rr_h_
#define _libint2_src_bin_libint_rr_h_

#include <string>
#include <vector>

using namespace std;


namespace libint2 {

  const int Libint2_DefaultVectorLength = 64;

  class GaussianShell {
    int l_;
    
    public:
    GaussianShell(int l) {l_ = l;};
    ~GaussianShell() {};
    
  };
  
  class CartesianGaussian : public GaussianShell {

    // cartesian exponents
    int n_[3];

    public:
    CartesianGaussian(int nx, int ny, int nz);
    ~CartesianGaussian() {};

  };

  class Operator {

    // Described name
    const std::string descr_;
    // short (<20 chars) ID label
    const std::string id_;
    
    // number of particles > 0
    char np_;
    // symmetry W.R.T. permutation of each pair of particles
    // 1 -- symmetric, -1 -- antisymmetric, 0 -- nonsymmetric
    // stored as a lower triangle (diagonal not included)
    vector<char> psymm_;
    
    public:
    Operator(const std::string& descr, const std::string& id, char np,
             const vector<char>& psymm);
    ~Operator();

    /// Returns full description of the operator
    const std::string& descr() const;
    /// Returns short label for the operator
    const std::string& id() const;
    
    /** Returns 1, 0, or -1, if the operator is symmetric, nonsymmetric,
        or antisymmetric with respect to permutation of particles i and j */
    const int psymm(int i, int j) const;

  };
  
  extern Operator TwoERep;
  extern Operator TwoEDist;

  class IntegralType {

    // Operator
    Operator O_;
    
    // number of basis functions in bra for each particle >= 0
    vector<char> nbra_;
    // number of basis functions in ket for each particle >= 0
    vector<char> nket_;
    
    public:
    IntegralType(const Operator& O, const vector<char>& nbra, const vector<char>& nket);
    ~IntegralType();
    
  };


  class IntegralSetType {

    public:
    IntegralSetType();
    virtual ~IntegralSetType();
    
    /// Number of elements (>0 if static, 0 if determined at runtime)
    virtual int num_elems();

  };

  class SingleIntegral : public IntegralSetType {

    public:
    SingleIntegral();

  };

  class IntegralShell : public IntegralSetType {

    public:
    IntegralShell();

  };

  class IntegralShellVector : public IntegralSetType {

    public:
    IntegralShellVector();

  };

  class IntegralVectorShell : public IntegralSetType {

    public:
    IntegralVectorShell();

  };

  class RecurrenceRelation {

    public:
    RecurrenceRelation();
    ~RecurrenceRelation();
    
    virtual const std::string cpp_function_name() =0;
    virtual const std::string cpp_source_name() =0;
    virtual const std::string cpp_header_name() =0;
    virtual std::ostream& cpp_source(std::ostream&) =0;

  };

  class DataTypedRR : public RecurrenceRelation {

    IntegralSetType ST_;

    public:
    DataTypedRR(const IntegralSetType&);
    ~DataTypedRR();

  };


  template<typename Q> class QuantumSet {

    Q quanta_;

public:
    QuantumSet(Q quanta);
    ~QuantumSet();

    /// Increment quantum number i
    void inc_quantum(unsigned int i) { quanta_.inc(i); };
    /// Decrement quantum number i
    void dec_quantum(unsigned int i) { quanta_.dec(i); };

  };

  template<typename T, unsigned int N> class QuantumNumbers {
    T qn_[N];

public:
    QuantumNumbers(T qn[N]);
    ~QuantumNumbers();

    /// Increment quantum number i
    void inc(unsigned int i) { ++qn[i]; };
    /// Decrement quantum number i
    void dec(unsigned int i) {
      if (qn[i] == T(0))
        throw std::runtime_error("QuantumNumber::dec -- quantum number already zero");
      --qn[i];
    };
  };

  template<typename Q>
    QuantumSet<Q>::QuantumSet(Q quanta) :
    quanta_(quanta)
  {
  };

  template<typename Q>
    QuantumSet<Q>::~QuantumSet()
  {
  };

  template<typename T, unsigned int N>
    QuantumNumbers<T,N>::QuantumNumbers(T qn[N])
  {
      for(int i=0; i<N; i++) qn_[i] = qn[i];
  };

  template<typename T, unsigned int N>
    QuantumNumbers<T,N>::~QuantumNumbers()
  {
  };

  typedef QuantumNumbers<unsigned int, 1> ShellQuantumNumbers;
  typedef QuantumNumbers<unsigned int, 3> CartGaussQuantumNumbers;
  typedef QuantumSet< ShellQuantumNumbers > GaussShell;
  typedef QuantumSet< CartGaussQuantumNumbers > CartGauss;
  
  
};

#endif

