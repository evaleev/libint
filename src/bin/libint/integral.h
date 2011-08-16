
#ifndef _libint2_src_bin_libint_integral_h_
#define _libint2_src_bin_libint_integral_h_

#include <smart_ptr.h>
#include <dgvertex.h>
#include <oper.h>
#include <iter.h>
#include <policy_spec.h>
#include <quanta.h>
#include <equiv.h>
#include <singl_stack.h>
#include <class_registry.h>
#include <global_macros.h>

#if USE_BRAKET_H
#  include <braket.h>
#endif

using namespace std;

extern long living_count;

namespace libint2 {

  /**
     This is an abstract base for sets of all types of integrals. Functions can be of any type
     derived from BasisFunctionSet.
  */
  template <class BasisFunctionSet> class IntegralSet {

  public:
    virtual ~IntegralSet() {};

#if USE_BRAKET_H
    /// Return the number of particles
    virtual unsigned int num_part() const =0;
    /// Return the number of functions for particle p
    virtual unsigned int num_func_bra(unsigned int p) const =0;
    /// Return the number of functions for particle p
    virtual unsigned int num_func_ket(unsigned int p) const =0;
    /// Obtain pointers to ith BasisFunctionSet for particle p in bra
    virtual const BasisFunctionSet& bra(unsigned int p, unsigned int i) const =0;
    /// Obtain pointers to ith BasisFunctionSet for particle p in ket
    virtual const BasisFunctionSet& ket(unsigned int p, unsigned int i) const =0;
    /// Obtain pointers to ith BasisFunctionSet for particle p in bra
    virtual BasisFunctionSet& bra(unsigned int p, unsigned int i) =0;
    /// Obtain pointers to ith BasisFunctionSet for particle p in ket
    virtual BasisFunctionSet& ket(unsigned int p, unsigned int i) =0;
#else
    /// Return the number of particles
    virtual unsigned int np() const =0;
    /// Obtain pointers to ith BasisFunctionSet for particle p in bra
    virtual const SafePtr<BasisFunctionSet> bra(unsigned int p, unsigned int i) const =0;
    /// Obtain pointers to ith BasisFunctionSet for particle p in ket
    virtual const SafePtr<BasisFunctionSet> ket(unsigned int p, unsigned int i) const =0;
#endif
  };

  /**
     GenIntegralSet is a set of integrals over functions derived from BFS.

     @tparam Oper An operator or a set of operators. Oper must be derived from OperSet.
     @tparam BraSetType Type describing a set of bra functions. An example of a class that can be used as BraSetType
     is VectorBraket.
     @tparam KetSetType Type describing a set of ket functions. An example of a class that can be used as KetSetType
     is VectorBraket.
     @tparam AuxQuanta Type describing a set of auxiliary quantum numbers. AuxQuanta should be derived from QuantumSet.
  */
  template <class Oper, class BFS, class BraSetType, class KetSetType, class AuxQuanta = EmptySet>
    class GenIntegralSet :
    public IntegralSet<BFS>, public DGVertex,
    public EnableSafePtrFromThis< GenIntegralSet<Oper,BFS,BraSetType,KetSetType,AuxQuanta> >
    {
      public:
      typedef GenIntegralSet this_type;
      /// GenIntegralSet is a set of these subobjects
      typedef GenIntegralSet<typename Oper::iter_type, BFS, typename BraSetType::iter_type, typename KetSetType::iter_type, typename AuxQuanta::iter_type> iter_type;
      /// GenIntegralSet is derived from IntegralSet
      typedef IntegralSet<BFS> parent_type;
      /// This type provides comparison operations on pointers to GenIntegralSet
      typedef PtrEquiv<GenIntegralSet> PtrComp;
#if  USE_INT_KEY_TO_HASH
      typedef mpz_class key_type;
#else
      typedef std::string key_type;
#endif
      /// This the type of the object that manages GenIntegralSet's as Singletons
      typedef SingletonStack<GenIntegralSet,key_type> SingletonManagerType;
      /// This is the type of the operator
      typedef Oper OperType;
      /// This is the real type of basis functions
      typedef typename BraSetType::bfs_type BasisFunctionType;

      /** No constructors are public since this is a singleton-like quantity.
          Instead, access is provided through Instance().

          Derived classes do not have to be Singletons -- hence protected constructors are provided.
      */

      /// GenIntegralSet objects can be handled through the base pointer -- hence the destructor is virtual
      virtual ~GenIntegralSet();

      /// Returns a pointer to a unique instance, a la Singleton
      static const SafePtr<GenIntegralSet> Instance(const BraSetType& bra, const KetSetType& ket, const AuxQuanta& aux, const Oper& oper=Oper());

      /// Comparison operator
      virtual bool operator==(const GenIntegralSet&) const;
      /// Specialization of DGVertex::equiv()
      bool equiv(const SafePtr<DGVertex>& v) const
      {
        return PtrComp::equiv(this,v);
      }
      /// Specialization of DGVertex::size()
      virtual const unsigned int size() const;
      /// Specialization of DGVertex::label()
      virtual const std::string& label() const;
      /// Specialization of DGVertex::id()
      virtual const std::string& id() const;
      /// Specialization of DGVertex::description()
      virtual std::string description() const;

      /// Implementation of IntegralSet::num_part
      unsigned int num_part() const { return OperType::Properties::np; }
      /// Implementation of IntegralSet::num_func_bra
      virtual unsigned int num_func_bra(unsigned int p) const { return bra_.num_members(p); }
      /// Implementation of IntegralSet::num_func_ket
      virtual unsigned int num_func_ket(unsigned int p) const { return ket_.num_members(p); }
      /// Implementation of IntegralSet::bra() const
      typename BraSetType::bfs_cref bra(unsigned int p, unsigned int i) const;
      /// Implementation of IntegralSet::ket() const
      typename KetSetType::bfs_cref ket(unsigned int p, unsigned int i) const;
      /// Implementation of IntegralSet::bra()
      typename BraSetType::bfs_ref bra(unsigned int p, unsigned int i);
      /// Implementation of IntegralSet::ket()
      typename KetSetType::bfs_ref ket(unsigned int p, unsigned int i);

      typedef BraSetType BraType;
      typedef KetSetType KetType;
      typedef Oper OperatorType;
      typedef AuxQuanta AuxQuantaType;

      /// Obtain the operator
      const SafePtr<Oper> oper() const;
      /// Obtain const ref to bra
      const BraType& bra() const;
      /// Obtain const ref to bra
      const KetType& ket() const;
      /// Obtain the auxiliary quanta
      const SafePtr<AuxQuanta> aux() const;

      /// Implements Hashable::key()
      key_type key() const {
        if (key_ == 0) compute_key();
        return key_;
      }

      /// Reimplements DGVertex::unregister()
      void unregister() const;

      protected:
      // Basic Integral constructor. It is protected so that derived classes don't have to behave like singletons
      GenIntegralSet(const Oper& oper, const BraSetType& bra, const KetSetType& ket, const AuxQuanta& aux);
      /// computes a key. it's protected so that derived classes can use it to implement smart caching in constructors
      static key_type compute_key(const Oper& O, const BraType& bra, const KetType& ket, const AuxQuanta& aux) {
#define TEST_KEYTYPE_SAFETY 0
#if TEST_KEYTYPE_SAFETY
        key_type remainder = UINT64_MAX;
        remainder /= (key_type)aux.max_key(); assert(remainder != 0);
        remainder /= (key_type)ket.max_key(); assert(remainder != 0);
        remainder /= (key_type)bra.max_key(); assert(remainder != 0);
        remainder /= (key_type)O.max_key; assert(remainder != 0);
#endif
        key_type key = ( (key_type(O.key()) * key_type(bra.max_key()) + key_type(bra.key()) ) * key_type(ket.max_key()) +
			 key_type(ket.key()) ) * key_type(aux.max_key()) + key_type(aux.key());
        return key;
      }

      BraSetType bra_;
      KetSetType ket_;

      /// set size to sz
      void set_size(unsigned int sz);
      /// Specialization of DGVertex::this_precomputed()
      virtual bool this_precomputed() const { return false; }
      /// Resets all cached values
      void reset_cache() { key_ = 0; size_ = 0; }

    private:
      //
      // All integrals are Singletons by nature, therefore they must be treated as such
      // 1) No public constructors are provided
      // 2) protected members are provided to implement Singleton-type functionality
      //
      GenIntegralSet();
      GenIntegralSet(const GenIntegralSet&);
      // Copy is not permitted
      GenIntegralSet& operator=(const GenIntegralSet& source);

      // This is used to manage GenIntegralSet objects as singletons
      static SingletonManagerType singl_manager_;

      // The operator needs to be a real object rather than real type to be able to construct a SubIterator, etc.
      SafePtr<Oper> O_;
      // Same for AuxQuanta
      SafePtr<AuxQuanta> aux_;

      // size of the integral set
      mutable unsigned int size_;

      // unique label -- no 2 GenIntegralSet instances can have the same label!
      mutable std::string label_;
      // generates label_
      std::string generate_label() const;
      // key
      mutable key_type key_;

      /// computes and caches key
      void compute_key() const {
        key_ = compute_key(*(O_.get()),bra_,ket_,*(aux_.get()));
      }

    };

#if USE_INT_KEY_TO_HASH
  // I use key() to hash GenIntegralSet. Therefore compute_key() must return unique keys!
  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    typename GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::SingletonManagerType
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::singl_manager_(&GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::key);
#else
  // I use label() to hash GenIntegralSet. Therefore labels must be unique!
  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    typename GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::SingletonManagerType
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::singl_manager_(&GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::label);
#endif

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::GenIntegralSet(const Op& oper, const BraSetType& bra, const KetSetType& ket, const AuxQuanta& aux) :
    DGVertex(ClassInfo<GenIntegralSet>::Instance().id()), bra_(bra), ket_(ket), O_(SafePtr<Op>(new Op(oper))), aux_(SafePtr<AuxQuanta>(new AuxQuanta(aux))),
    size_(0), label_()
    {
      if (Op::Properties::np != bra.num_part())
        throw std::runtime_error("GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::GenIntegralSet(bra,ket) -- number of particles in bra doesn't match that in the operator");
      if (Op::Properties::np != ket.num_part())
        throw std::runtime_error("GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::GenIntegralSet(bra,ket) -- number of particles in ket doesn't match that in the operator");
      compute_key();
#if DEBUG
      std::cout << "GenIntegralSet: constructed " << label() << std::endl;
      std::cout << "GenIntegralSet: living_count = " << ++living_count << std::endl;
#endif
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::~GenIntegralSet()
    {
#if DEBUG
      std::cout << "GenIntegralSet: destructed " << label() << std::endl;
      std::cout << "GenIntegralSet: living_count = " << --living_count << std::endl;
#endif
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const SafePtr< GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta> >
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::Instance(const BraSetType& bra, const KetSetType& ket, const AuxQuanta& aux, const Op& oper)
    {
      typedef typename SingletonManagerType::value_type map_value_type;
      key_type key = compute_key(oper,bra,ket,aux);
      const map_value_type& val = singl_manager_.find(key);
      if (!val.second) {
	SafePtr<this_type> this_int(new this_type(oper,bra,ket,aux));
	// Use singl_manager_ to make sure this is a new object of this type
	const map_value_type& val = singl_manager_.find(this_int);
	val.second->instid_ = val.first;
	this_int.reset();
	return val.second;
      }
      return val.second;
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    typename BraSetType::bfs_cref
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::bra(unsigned int p, unsigned int i) const
    {
      return bra_.member(p,i);
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    typename KetSetType::bfs_cref
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::ket(unsigned int p, unsigned int i) const
    {
      return ket_.member(p,i);
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    typename BraSetType::bfs_ref
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::bra(unsigned int p, unsigned int i)
    {
      reset_cache();
      return bra_.member(p,i);
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    typename KetSetType::bfs_ref
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::ket(unsigned int p, unsigned int i)
    {
      reset_cache();
      return ket_.member(p,i);
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const typename GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::BraType&
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::bra() const
    {
      return bra_;
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const typename GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::KetType&
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::ket() const
    {
      return ket_;
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const SafePtr<Op>
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::oper() const
    {
      return O_;
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const SafePtr<AuxQuanta>
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::aux() const
    {
      return aux_;
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    bool
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::operator==(const this_type& a) const
    {
#if USE_INT_KEY_TO_COMPARE
      return key() == a.key();
#else
      bool aux_equiv = PtrEquiv<AuxQuanta>::equiv(aux_,a.aux_);
      if (!aux_equiv) return false;
      bool oper_equiv = PtrEquiv<Op>::equiv(O_,a.O_);
      bool bra_equiv = PtrEquiv<BraSetType>::equiv(bra_,a.bra_);
      if (!bra_equiv) return false;
      bool ket_equiv = PtrEquiv<KetSetType>::equiv(ket_,a.ket_);
      if (!ket_equiv) return false;
      return true;
#endif
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    void
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::unregister() const
    {
      singl_manager_.remove(const_pointer_cast<this_type,const this_type>(EnableSafePtrFromThis<this_type>::SafePtr_from_this()));
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const unsigned int
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::size() const
    {
      if (size_ > 0)
        return size_;

#if COMPUTE_SIZE_DIRECTLY
    // WARNING: will not work if aux_ includes "sets" of numbers, i.e. when a quantum number can be a set of numbers
    // but I don't at the moment see why this would be needed
    size_ = bra_.size() * ket_.size() * O_->num_oper();
#else
      // compute size
      SafePtr<this_type> this_ptr = const_pointer_cast<this_type,const this_type>(EnableSafePtrFromThis<GenIntegralSet>::SafePtr_from_this());
      SafePtr< SubIteratorBase<this_type> > siter(new SubIteratorBase<this_type>(this_ptr));
      size_ = siter->num_iter();
      if (size_ == 0)
        throw std::runtime_error("GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::size() -- size is 0");
#endif

      return size_;
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    void
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::set_size(unsigned int sz)
    {
      size_ = sz;
    }

  template <class BraSetType, class KetSetType, class AuxQuanta, class Op>
    std::string
    genintegralset_label(const BraSetType& bra, const KetSetType& ket, const SafePtr<AuxQuanta>& aux, const SafePtr<Op>& O)
    {
      ostringstream os;
      os << "< ";
      for(unsigned int p=0; p<Op::Properties::np; p++) {
        unsigned int nbra = bra.num_members(p);
        for(unsigned int i=0; i<nbra; i++)
#if USE_BRAKET_H
          os << bra.member(p,i).label() << "(" << p << ") ";
#else
          os << bra.member(p,i)->label() << "(" << p << ") ";
#endif
      }
      os << "| " << O->label() << " | ";
      for(unsigned int p=0; p<Op::Properties::np; p++) {
        unsigned int nket = ket.num_members(p);
        for(unsigned int i=0; i<nket; i++)
#if USE_BRAKET_H
          os << ket.member(p,i).label() << "(" << p << ") ";
#else
          os << ket.member(p,i)->label() << "(" << p << ") ";
#endif
      }
      os << "> ^ { " << aux->label() << " }";
      return os.str();
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const std::string&
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::label() const
    {
      if (label_.empty())
        label_ = generate_label();
      return label_;
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    std::string
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::generate_label() const
    {
      return genintegralset_label(bra_,ket_,aux_,O_);
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    const std::string&
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::id() const
    {
      return label();
    }

  template <class Op, class BFS, class BraSetType, class KetSetType, class AuxQuanta>
    std::string
    GenIntegralSet<Op,BFS,BraSetType,KetSetType,AuxQuanta>::description() const
    {
    ostringstream os;
    os << " GenIntegralSet: " << label();
    const std::string descr = os.str();
    return descr;
    }

};

#endif
