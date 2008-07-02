
#ifndef _libint2_src_bin_libint_oper_h_
#define _libint2_src_bin_libint_oper_h_

#include <string>
#include <hashable.h>
#include <global_macros.h>
#include <util.h>

namespace libint2 {

  /** Permutational symmetries: antisymmetric(anti), symmetric(symm), nonsymmetric (nonsymm),
      some more complicated symmetry (nonstd) */
  typedef struct {
    typedef enum {anti=-1, symm=1, nonsymm=0, nonstd=-2} type;
  } PermutationalSymmetry;

  /** OperatorProperties describes various properties of an operator or operator set
      np -- number of particles
      multi -- true if multiplicative
      psymmetry -- symmetry with respect to permutation of bra and ket
    */
  template <unsigned int NP, bool multi, PermutationalSymmetry::type psymmetry>
    class OperatorProperties {
    public:
      static const unsigned int np = NP;
      static const bool multiplicative = multi;
      static const PermutationalSymmetry::type psymm = psymmetry;
    };

  /** OperSet is the base class for all (sets of) operators.
     OperSet's must be constructable using
     SafePtr<OperSet> or SafePtr<ConstructablePolymorphically>.
  */
  class OperSet : public ConstructablePolymorphically {
    public:
      virtual ~OperSet() {};

      /// Returns full human-readable description of the operator
      virtual std::string description() const =0;
      /// Returns short label for the operator
      virtual std::string label() const =0;
      /** Returns 1, 0, or -1, if each operator in the set is symmetric, nonsymmetric,
        or antisymmetric with respect to permutation of particles i and j */
      virtual int psymm(int i, int j) const =0;
      /** Returns 1, 0, or -1, if each operator in the set is Hermitian, non-Hermitian,
        or anti-Hermitian w.r.t. particle p */
      virtual int hermitian(int p) const =0;

      /// Number of operators in the set
      virtual const unsigned int num_oper() const =0;
    };

  /** Oper is OperSet characterized by properties Props
  */
  template <class Props>
    class Oper : public OperSet {
    public:
      typedef Props Properties;
      virtual ~Oper() {}

      /// Implementation of OperSet::psymm()
      int psymm(int i, int j) const;
      /// Implementation of OperSet::hermitian()
      int hermitian(int p) const;

      bool operator==(const Oper&) const;

    protected:
      /// The only declared constructor is only useable by derived classes
      Oper() {}

    private:
      /// Must be overloaded by derived classes if Props::psymm == PermutationalSymmetry::nonstd
      virtual int nonstd_psymm(int i, int j) const { throw ProgrammingError("nonstd_psymm is not overloaded"); }
      /// Must be overloaded by derived classes if Props::multi != true
      virtual int nonstd_hermitian(int p) const { throw ProgrammingError("nonstd_hermitian is not overloaded"); }
  };

  template <class Props>
    int
    Oper<Props>::psymm(int i, int j) const
    {
      if (i<0 || i>=Props::np)
        throw std::runtime_error("Oper<Props>::psymm(i,j) -- index i out of bounds");
      if (j<0 || j>=Props::np)
        throw std::runtime_error("Oper<Props>::psymm(i,j) -- index j out of bounds");
      if (i == j)
        return 1;

      switch(Props::psymm) {
        case PermutationalSymmetry::anti:
        return -1;
        case PermutationalSymmetry::symm:
        return 1;
        case PermutationalSymmetry::nonsymm:
        return 0;
        case PermutationalSymmetry::nonstd:
        return nonstd_psymm(i,j);
      }
    }

  template <class Props>
    int
    Oper<Props>::hermitian(int p) const
    {
      if (Props::multiplicative)
        return +1;
      else
        return nonstd_hermitian(p);
    }

  template <class Props>
  bool
    Oper<Props>::operator==(const Oper& a) const
    {
      return true;
    }

//////////////////////////////

  /** GenOper is a single operator described by descriptor Descr
  */
  template <class Descr>
    class GenOper : public Oper<typename Descr::Properties>, public Hashable<unsigned,ComputeKey> {
    public:
      typedef typename Descr::Properties Properties;
      typedef Oper<Properties> parent_type;
      /// GenOper is not a set
      typedef GenOper iter_type;

      const unsigned int num_oper() const { return 1; };
      /// Implementation of Hashable::key()
      unsigned int key() const { return descr_.key(); }
      /// Range of key is [0,Descr::max_key)
      static const unsigned int max_key = Descr::max_key;
      /// Implementation of OperSet::description()
      std::string description() const { return descr_.description(); }
      /// Implementation of OperSet::label()
      std::string label() const { return descr_.label(); }
      /// Return the descriptor object
      Descr descr() const { return descr_; }

      GenOper(Descr descr = Descr()) : descr_(descr) {}
      GenOper(const SafePtr<GenOper>& o) : descr_(o->descr_) {}
      GenOper(const SafePtr<OperSet>& o) : descr_(require_dynamic_cast<GenOper,OperSet>(o)->descr_) {}
      GenOper(const SafePtr<ConstructablePolymorphically>& o) : descr_(require_dynamic_cast<GenOper,ConstructablePolymorphically>(o)->descr_) {}
      GenOper(const ConstructablePolymorphically& o) : descr_(require_dynamic_cast<GenOper,ConstructablePolymorphically>(&o)->descr_) {}
      ~GenOper() {}

    private:
      Descr descr_;

      /// Implementation of Oper::nonstd_psymm()
      int nonstd_psymm(int i, int j) const {
        // TODO: figure out how to call this only if Desc::Properties::psymm == PermutationalSymmetry::nonstd
        if (Descr::Properties::psymm == PermutationalSymmetry::nonstd)
          return descr_.psymm(i,j);
        else throw ProgrammingError("GenOper::nonstd_psymm -- descriptor is not nonstd");
      }

      /// Implementation of Oper::nonstd_hermitian()
      int nonstd_hermitian(int i) const {
        // TODO: figure out how to call this only if Desc::Properties::psymm == PermutationalSymmetry::nonstd
        if (!Descr::Properties::multiplicative)
          return descr_.hermitian(i);
        else throw ProgrammingError("GenOper::nonstd_hermitian -- this operator is multiplicative");
      }

    };

//////////////////////////////

  typedef OperatorProperties<2,true,PermutationalSymmetry::symm> MultiplicativeSymm2Body_Props;
  typedef OperatorProperties<2,true,PermutationalSymmetry::nonsymm> MultiplicativeNonsymm2Body_Props;
  typedef OperatorProperties<2,false,PermutationalSymmetry::symm> NonmultiplicativeSymm2Body_Props;
  typedef OperatorProperties<2,false,PermutationalSymmetry::nonsymm> NonmultiplicativeNonsymm2Body_Props;

  /** GenMultSymmOper is a generic multiplicative symmetric N-body operator
  */
  template <unsigned int N>
  struct GenMultSymmOper_Descr {
    typedef OperatorProperties<N,true,PermutationalSymmetry::symm> Properties;
    static const unsigned int max_key = 1;
    unsigned int key() const { return 0; }
    std::string description() const { return "generic multiplicative symmetric operator"; }
    std::string label() const { return "GenMultSymmOper"; }
    int psymm(int i, int j) const;
    int hermitian(int i) const;
  };
  typedef GenOper< GenMultSymmOper_Descr<2>  > GenMultSymm2BodyOper;

  /** TwoPRep is the two-body repulsion operator.
  */
  struct TwoPRep_Descr {
    typedef MultiplicativeSymm2Body_Props Properties;
    static const unsigned int max_key = 1;
    unsigned int key() const { return 0; }
    std::string description() const { return "1/r_{12}"; }
    std::string label() const { return "TwoPRep"; }
    int psymm(int i, int j) const;
    int hermitian(int i) const;
  };
  typedef GenOper<TwoPRep_Descr> TwoPRep;

  /** R12_k_G12 is a two-body operator of form r_{12}^k * exp(-\gamma * r_{12}),
      where k is an integer and \gamma is a positive real number.
  */
  class R12_k_G12_Descr {
  public:
    typedef MultiplicativeSymm2Body_Props Properties;
    /// K can range from -1 to 4
    R12_k_G12_Descr(int K) : K_(K) { assert(K >= -1 && K <= 4); }
    R12_k_G12_Descr(const R12_k_G12_Descr& a) : K_(a.K_) {}
    static const unsigned int max_key = 5;
    unsigned int key() const { return K_ + 1; }
    std::string description() const { return label_(K_); }
    std::string label() const { return symbol_(K_); }
    int K() const { return K_; }
    int psymm(int i, int j) const;
    int hermitian(int i) const;
  private:
    R12_k_G12_Descr();
    static std::string label_(int K);
    static std::string symbol_(int K);
    int K_;
  };
  typedef GenOper<R12_k_G12_Descr> R12kG12;

  /** Ti_G12 is a two-body operator of form [T_i, G12],
      where i is particle index (0 or 1) and G12 is a Gaussian Geminal.
  */
  class Ti_G12_Descr {
  public:
    typedef NonmultiplicativeNonsymm2Body_Props Properties;
    /// K can range from 0 to 1
    static const unsigned int max_key = 2;
    Ti_G12_Descr(int K) : K_(K) { assert(K >= 0 && K <= 1); }
    Ti_G12_Descr(const Ti_G12_Descr& a) : K_(a.K_) {}
    unsigned int key() const { return K_; }
    std::string description() const { return label_(K_); }
    std::string label() const { return symbol_(K_); }
    int K() const { return K_; }
    int psymm(int i, int j) const;
    int hermitian(int i) const { if (i != K_) return +1; else return -1; }
  private:
    Ti_G12_Descr();
    static std::string label_(int K);
    static std::string symbol_(int K);
    int K_;
  };
  typedef GenOper<Ti_G12_Descr> TiG12;

  /** r_1.r_1 x g12 is a result of differentiation of exp( - a r_1^2 - a r_2^2 - c r_{12}^2) geminal .
  */
  struct R1dotR1_G12_Descr {
    typedef MultiplicativeNonsymm2Body_Props Properties;
    static const unsigned int max_key = 1;
    unsigned int key() const { return 0; }
    std::string description() const { return "r_1.r_1 x G12"; }
    std::string label() const { return "R1dotR1_G12"; }
    int psymm(int i, int j) const;
    int hermitian(int i) const;
  };
  typedef GenOper< R1dotR1_G12_Descr > R1dotR1_G12;

  /** r_2.r_2 x g12 is a result of differentiation of exp( - a r_1^2 - a r_2^2 - c r_{12}^2) geminal .
  */
  struct R2dotR2_G12_Descr {
    typedef MultiplicativeNonsymm2Body_Props Properties;
    static const unsigned int max_key = 1;
    unsigned int key() const { return 0; }
    std::string description() const { return "r_2.r_2 x G12"; }
    std::string label() const { return "R2dotR2_G12"; }
    int psymm(int i, int j) const;
    int hermitian(int i) const;
  };
  typedef GenOper< R2dotR2_G12_Descr > R2dotR2_G12;

  /** r_1.r_2 x g12 is a result of differentiation of exp( - a r_1^2 - a r_2^2 - c r_{12}^2) geminal .
  */
  struct R1dotR2_G12_Descr {
    typedef MultiplicativeSymm2Body_Props Properties;
    static const unsigned int max_key = 1;
    unsigned int key() const { return 0; }
    std::string description() const { return "r_1.r_2 x G12"; }
    std::string label() const { return "R1dotR2_G12"; }
    int psymm(int i, int j) const;
    int hermitian(int i) const;
  };
  typedef GenOper< R1dotR2_G12_Descr > R1dotR2_G12;

  /** \f$ (\nabla_1 \cdot g_{12}') (g_{12}' \cdot \nabla_1) \f$ is a component of \f$ [g_{12}, [\nabla_1^4, g_{12}]] \f$ integral .
  */
  struct DivG12prime_xTx_Descr {
    typedef NonmultiplicativeNonsymm2Body_Props Properties;
    static const unsigned int max_key = 1;
    unsigned int key() const { return 0; }
    std::string description() const { return "(\\nabla \\cdot g_{12}') (g_{12}' \\cdot \\nabla)"; }
    std::string label() const { return "DivG12prime_xTx"; }
    int psymm(int i, int j) const;
    int hermitian(int i) const { if (i != 0) return +1; else return -1; }
  };
  typedef GenOper< DivG12prime_xTx_Descr > DivG12prime_xTx;

};

#endif

