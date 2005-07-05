
#include <integral.h>

#ifndef _libint2_src_bin_libint_r12kg121111_h_
#define _libint2_src_bin_libint_r12kg121111_h_

using namespace std;

namespace libint2 {
  
  /** R12kG12_11_11_base is the base for all 2-body integral over the R12_k_G12
      operator with one basis function for each particle in bra and ket
    */
  template <class BFSet>
    struct R12kG12_11_11_base {
      virtual bool instance_of_R12kG12_11_11() =0;
    };

  /**
     mType is the type that describes the auxiliary index of
     any 2-body integral.
  */
  typedef QuantumNumbers<unsigned int,1> mType;
  /**
     Most basic type -- R12kG12_11_11 --
     has one bfs for each particle in bra and ket.
     Note that GenIntegralSet is initialized with an abstract type libint2::BFSet,
     from which BFS derives.
  */
  template <class BFS, int K> class R12kG12_11_11 :
    public GenIntegralSet< R12_k_G12<K>, IncableBFSet, VectorBraket<BFS>, VectorBraket<BFS>, mType >,
    public R12kG12_11_11_base<BFS>
    {
    public:
      typedef R12_k_G12<K> OperType;
      typedef VectorBraket<BFS> BraType;
      typedef VectorBraket<BFS> KetType;
      typedef mType AuxIndexType;
      typedef R12kG12_11_11 this_type;
      /// This the type of the object that manages GenIntegralSet's as Singletons
      typedef SingletonStack<R12kG12_11_11,std::string> SingletonManagerType;
      
      /// R12kG12_11_11 is a set of these subobjects
      typedef R12kG12_11_11<typename BFS::iter_type,K> iter_type;
      /// This is the immediate parent
      typedef GenIntegralSet< OperType, IncableBFSet, VectorBraket<BFS>, VectorBraket<BFS>, AuxIndexType > parent_type;
      /// This class provides comparison operations on pointers
      typedef PtrEquiv<this_type> PtrComp;

      /* This "constructor" takes basis function sets, in Mulliken ordering.
         Returns a pointer to a unique instance, a la Singleton
      */
      static const SafePtr<R12kG12_11_11> Instance(const BFS& bra0, const BFS& ket0, const BFS& bra1, const BFS& ket1, unsigned int m);
      /// Returns a pointer to a unique instance, a la Singleton
      static const SafePtr<R12kG12_11_11> Instance(const VectorBraket<BFS>& bra, const VectorBraket<BFS>& ket, const AuxIndexType& aux);
      
      unsigned int m() const { return parent_type::aux()->elem(0); };

      /// Comparison operator
      bool operator==(const this_type&) const;
#if OVERLOAD_GENINTEGRALSET_LABEL
      /// Specialization of GenIntegralSet::label()
      const std::string& label() const;
#endif

    private:
      // This constructor is also private and not implemented since all Integral's are Singletons. Use Instance instead.
      R12kG12_11_11(const VectorBraket<BFS>& bra, const VectorBraket<BFS>& ket, const AuxIndexType& aux);

      // This is used to manage GenIntegralSet objects as singletons
      static SingletonManagerType singl_manager_;

      /// Implements DGVertex::this_precomputed()
      bool this_precomputed() const;
      /// Implements R12kG12_11_11_base::instance_of_R12kG12_11_11()
      bool instance_of_R12kG12_11_11() { return true; }
#if OVERLOAD_GENINTEGRALSET_LABEL
      mutable std::string label_;
#endif
    };

  // I use label() to hash R12kG12_11_11. Therefore labels must be unique!
  template <class BFS, int K>
    typename R12kG12_11_11<BFS,K>::SingletonManagerType
    R12kG12_11_11<BFS,K>::singl_manager_(&R12kG12_11_11<BFS,K>::label);
  
  template <class BFS, int K>
    R12kG12_11_11<BFS,K>::R12kG12_11_11(const VectorBraket<BFS>& bra, const VectorBraket<BFS>& ket,  const AuxIndexType& aux) :
    parent_type(R12_k_G12<K>(),bra, ket, aux)
    {
      if (bra.num_members(0) != 1)
        throw std::runtime_error("R12kG12_11_11<BFS,K>::R12kG12_11_11(bra,ket) -- number of BFSs in bra for particle 0 must be 1");
      if (bra.num_members(1) != 1)
        throw std::runtime_error("R12kG12_11_11<BFS,K>::R12kG12_11_11(bra,ket) -- number of BFSs in bra for particle 1 must be 1");
      if (ket.num_members(0) != 1)
        throw std::runtime_error("R12kG12_11_11<BFS,K>::R12kG12_11_11(bra,ket) -- number of BFSs in ket for particle 0 must be 1");
      if (ket.num_members(1) != 1)
        throw std::runtime_error("R12kG12_11_11<BFS,K>::R12kG12_11_11(bra,ket) -- number of BFSs in ket for particle 1 must be 1");
    }

  template <class BFS, int K>
    const SafePtr< R12kG12_11_11<BFS,K> >
    R12kG12_11_11<BFS,K>::Instance(const VectorBraket<BFS>& bra, const VectorBraket<BFS>& ket, const AuxIndexType& aux)
    {
      SafePtr<R12kG12_11_11> this_int(new R12kG12_11_11<BFS,K>(bra,ket,aux));
      // Use singl_manager_ to make sure this is a new object of this type
      const typename SingletonManagerType::value_type& val = singl_manager_.find(this_int);
      val.second->instid_ = val.first;
      return val.second;
    }

  template <class BFS, int K>
    const SafePtr< R12kG12_11_11<BFS,K> >
    R12kG12_11_11<BFS,K>::Instance(const BFS& bra0, const BFS& ket0, const BFS& bra1, const BFS& ket1, unsigned int m)
    {
#if USE_BRAKET_H
      typedef BFS BFSRef;
      BFSRef bra0_ref(bra0);
      BFSRef bra1_ref(bra1);
      BFSRef ket0_ref(ket0);
      BFSRef ket1_ref(ket1);
#else
      typedef SafePtr<BFS> BFSRef;
      BFSRef bra0_ref(new BFS(bra0));
      BFSRef bra1_ref(new BFS(bra1));
      BFSRef ket0_ref(new BFS(ket0));
      BFSRef ket1_ref(new BFS(ket1));
#endif
      vector<BFSRef> vbra0;  vbra0.push_back(bra0_ref);
      vector<BFSRef> vbra1;  vbra1.push_back(bra1_ref);
      vector<BFSRef> vket0;  vket0.push_back(ket0_ref);
      vector<BFSRef> vket1;  vket1.push_back(ket1_ref);
      vector< vector<BFSRef> > vvbra;  vvbra.push_back(vbra0);  vvbra.push_back(vbra1);
      vector< vector<BFSRef> > vvket;  vvket.push_back(vket0);  vvket.push_back(vket1);
      VectorBraket<BFS> bra(vvbra);
      VectorBraket<BFS> ket(vvket);
      AuxIndexType aux(vector<unsigned int>(1,m));
      return Instance(bra,ket,aux);
    }
  
  template <class BFS, int K>
    bool
    R12kG12_11_11<BFS,K>::operator==(const R12kG12_11_11<BFS,K>& a) const
    {
      return parent_type::PtrComp::equiv(static_cast<const parent_type*>(this),a);
    }

#if OVERLOAD_GENINTEGRALSET_LABEL
  template <class BFS, int K>
    const std::string&
    R12kG12_11_11<BFS,K>::label() const
    {
      if (label_.empty()) {
	ostringstream os;
	os << "(" << parent_type::bra_.member(0,0)->label() << " "
	   << parent_type::ket_.member(0,0)->label()
	   << " | r_{12}^" << K << " * G12 | "
	   << parent_type::bra_.member(1,0)->label() << " "
	   << parent_type::ket_.member(1,0)->label() << ")^{" << m() <<"}";
	label_ = os.str();
      }
      return label_;
    };
#endif
  
  template <class BFS, int K>
    bool
    R12kG12_11_11<BFS,K>::this_precomputed() const
    {
      if (TrivialBFSet<BFS>::result == false)
        return false;
      else {
#if USE_BRAKET_H
        if (parent_type::bra_.member(0,0).zero() && parent_type::bra_.member(1,0).zero() &&
            parent_type::ket_.member(0,0).zero() && parent_type::ket_.member(1,0).zero())
#else
        if (parent_type::bra_.member(0,0)->zero() && parent_type::bra_.member(1,0)->zero() &&
            parent_type::ket_.member(0,0)->zero() && parent_type::ket_.member(1,0)->zero())
#endif
          return true;
        else
         return false;
      }
    }
  
  /** This is a collection of utility functions to extract
      compile-time information about instances of R12kG12_11_11 at runtime.
    */
  namespace R12kG12_11_11_Util {

    template <class BFS>
      int k(const SafePtr< R12kG12_11_11_base<BFS> >& ints) {
        typedef R12kG12_11_11_base<BFS> base_type;
        // try k=2
        {
          typedef R12kG12_11_11<BFS,2> inst_type;
          SafePtr<inst_type> iptr = dynamic_pointer_cast<inst_type,base_type>(ints);
          if (iptr != 0)
            return 2;
        }
        // try k=0
        {
          typedef R12kG12_11_11<BFS,0> inst_type;
          SafePtr<inst_type> iptr = dynamic_pointer_cast<inst_type,base_type>(ints);
          if (iptr != 0)
            return 0;
        }
        // try k=-1
        {
          typedef R12kG12_11_11<BFS,-1> inst_type;
          SafePtr<inst_type> iptr = dynamic_pointer_cast<inst_type,base_type>(ints);
          if (iptr != 0)
            return -1;
        }
        throw std::logic_error("R12kG12_11_11_Util::k<BFS>() -- unknown value of k");
      };

  };

  /// R12_m1_G12_11_11_... is a set of integrals over g12/r12
  typedef R12kG12_11_11<CGShell,-1> R12_m1_G12_11_11_sq;
  typedef R12kG12_11_11<CGF,-1> R12_m1_G12_11_11_int;

  /// R12_0_G12_11_11_... is a set of integrals over g12
  typedef R12kG12_11_11<CGShell,0> R12_0_G12_11_11_sq;
  typedef R12kG12_11_11<CGF,0> R12_0_G12_11_11_int;

};

#endif

