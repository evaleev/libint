/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint compiler.
 *
 *  Libint compiler is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint compiler is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint compiler.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_bin_libint_r2dotr2g121111_h_
#define _libint2_src_bin_libint_r2dotr2g121111_h_

#include <integral.h>

namespace libint2 {

#if 0
  /**
     R2dotR2G12_11_11 --
     integral over R2dotR2_G12 operator with one bfs for each particle in bra and ket.
  */
  template <class BFS> class R2dotR2G12_11_11 :
    public GenIntegralSet< R2dotR2_G12, IncableBFSet, typename DefaultTwoPBraket<BFS>::Result, typename DefaultTwoPBraket<BFS>::Result, EmptySet >
    {
    public:
      typedef BFS BasisFunctionType;
      typedef R2dotR2_G12 OperType;
      typedef typename DefaultTwoPBraket<BFS>::Result BraType;
      typedef typename DefaultTwoPBraket<BFS>::Result KetType;
      typedef EmptySet AuxIndexType;
      typedef R2dotR2G12_11_11<BFS> this_type;

      /// R2dotR2G12_11_11 is a set of these subobjects
      typedef R2dotR2G12_11_11<typename BFS::iter_type> iter_type;
      /// This is the immediate parent
      typedef GenIntegralSet< OperType, IncableBFSet, BraType, KetType, AuxIndexType > parent_type;
      /// This class provides comparison operations on pointers
      typedef PtrEquiv<this_type> PtrComp;

      typedef typename parent_type::key_type key_type;
      /// This the type of the object that manages GenIntegralSet's as Singletons
      typedef SingletonStack<R2dotR2G12_11_11,key_type> SingletonManagerType;

      /* This "constructor" takes basis function sets, in Mulliken ordering.
         Returns a pointer to a unique instance, a la Singleton
      */
      static const std::shared_ptr<R2dotR2G12_11_11> Instance(const BFS& bra0, const BFS& ket0, const BFS& bra1, const BFS& ket1);
      /// Returns a pointer to a unique instance, a la Singleton
      static const std::shared_ptr<R2dotR2G12_11_11> Instance(const BraType& bra, const KetType& ket, const AuxIndexType& aux);
      ~R2dotR2G12_11_11();

      /// Comparison operator
      bool operator==(const this_type&) const;
#if OVERLOAD_GENINTEGRALSET_LABEL
      /// Specialization of GenIntegralSet::label()
      const std::string& label() const;
#endif

    private:
      // This constructor is also private and not implemented since all Integral's are Singletons. Use Instance instead.
      R2dotR2G12_11_11(const BraType& bra, const KetType& ket, const AuxIndexType& aux);

      // This is used to manage GenIntegralSet objects as singletons
      static SingletonManagerType singl_manager_;

      /// Implements DGVertex::this_precomputed()
      bool this_precomputed() const;
#if OVERLOAD_GENINTEGRALSET_LABEL
      mutable std::string label_;
#endif

    };

#if USE_INT_KEY_TO_HASH
  template <class BFS>
    typename R2dotR2G12_11_11<BFS>::SingletonManagerType
    R2dotR2G12_11_11<BFS>::singl_manager_(&R2dotR2G12_11_11<BFS>::key);
#else
#error "USE_INT_KEY_TO_HASH must be set"
#endif

  template <class BFS>
    R2dotR2G12_11_11<BFS>::R2dotR2G12_11_11(const BraType& bra, const KetType& ket,  const AuxIndexType& aux) :
    parent_type(R2dotR2_G12(),bra, ket, aux)
    {
      if (bra.num_members(0) != 1)
        throw std::runtime_error("R2dotR2G12_11_11<BFS>::R2dotR2G12_11_11(bra,ket) -- number of BFSs in bra for particle 0 must be 1");
      if (bra.num_members(1) != 1)
        throw std::runtime_error("R2dotR2G12_11_11<BFS>::R2dotR2G12_11_11(bra,ket) -- number of BFSs in bra for particle 1 must be 1");
      if (ket.num_members(0) != 1)
        throw std::runtime_error("R2dotR2G12_11_11<BFS>::R2dotR2G12_11_11(bra,ket) -- number of BFSs in ket for particle 0 must be 1");
      if (ket.num_members(1) != 1)
        throw std::runtime_error("R2dotR2G12_11_11<BFS>::R2dotR2G12_11_11(bra,ket) -- number of BFSs in ket for particle 1 must be 1");
#if DEBUG
      std::cout << "R2dotR2G12_11_11: constructed " << this->label() << std::endl;
#endif
    }

  template <class BFS>
    R2dotR2G12_11_11<BFS>::~R2dotR2G12_11_11()
    {
#if DEBUG
      std::cout << "R2dotR2G12_11_11: destructed " << this->label() << std::endl;
#endif
    }

  template <class BFS>
    const std::shared_ptr< R2dotR2G12_11_11<BFS> >
    R2dotR2G12_11_11<BFS>::Instance(const BraType& bra, const KetType& ket, const AuxIndexType& aux)
    {
      typedef typename SingletonManagerType::value_type map_value_type;
      key_type key = compute_key(OperType(),bra,ket,aux);
      const map_value_type& val = singl_manager_.find(key);
      if (!val.second) {
std::shared_ptr<R2dotR2G12_11_11> this_int(new R2dotR2G12_11_11<BFS>(bra,ket,aux));
// Use singl_manager_ to make sure this is a new object of this type
const typename SingletonManagerType::value_type& val = singl_manager_.find(this_int);
val.second->instid_ = val.first;
return val.second;
      }
      return val.second;
    }

  template <class BFS>
    const std::shared_ptr< R2dotR2G12_11_11<BFS> >
    R2dotR2G12_11_11<BFS>::Instance(const BFS& bra0, const BFS& ket0, const BFS& bra1, const BFS& ket1)
    {
#if USE_BRAKET_H
      typedef BFS BFSRef;
      BFSRef bra0_ref(bra0);
      BFSRef bra1_ref(bra1);
      BFSRef ket0_ref(ket0);
      BFSRef ket1_ref(ket1);
#else
      typedef std::shared_ptr<BFS> BFSRef;
      BFSRef bra0_ref(new BFS(bra0));
      BFSRef bra1_ref(new BFS(bra1));
      BFSRef ket0_ref(new BFS(ket0));
      BFSRef ket1_ref(new BFS(ket1));
#endif
      std::vector<BFSRef> vbra0;  vbra0.push_back(bra0_ref);
      std::vector<BFSRef> vbra1;  vbra1.push_back(bra1_ref);
      std::vector<BFSRef> vket0;  vket0.push_back(ket0_ref);
      std::vector<BFSRef> vket1;  vket1.push_back(ket1_ref);
      std::vector< std::vector<BFSRef> > vvbra;  vvbra.push_back(vbra0);  vvbra.push_back(vbra1);
      std::vector< std::vector<BFSRef> > vvket;  vvket.push_back(vket0);  vvket.push_back(vket1);
      BraType bra(vvbra);
      KetType ket(vvket);
      AuxIndexType aux(std::vector<int>(0));
      return Instance(bra,ket,aux);
    }

  template <class BFS>
    bool
    R2dotR2G12_11_11<BFS>::operator==(const R2dotR2G12_11_11<BFS>& a) const
    {
      return parent_type::PtrComp::equiv(static_cast<const parent_type*>(this),a);
    }

#if OVERLOAD_GENINTEGRALSET_LABEL
  template <class BFS>
    const std::string&
    R2dotR2G12_11_11<BFS>::label() const
    {
      if (label_.empty()) {
ostringstream os;
os << "(" << parent_type::bra_.member(0,0)->label() << " "
   << parent_type::ket_.member(0,0)->label()
   << " | r_2^2 * G12 | "
   << parent_type::bra_.member(1,0)->label() << " "
   << parent_type::ket_.member(1,0)->label() << ")";
label_ = os.str();
      }
      return label_;
    };
#endif

  template <class BFS>
    bool
    R2dotR2G12_11_11<BFS>::this_precomputed() const
    {
      return false;
    }

  /// the following typedefs are useful
  typedef R2dotR2G12_11_11<CGShell> R2dotR2G12_11_11_sq;
  typedef R2dotR2G12_11_11<CGF> R2dotR2G12_11_11_int;
#endif
};  // namespace libint2

#endif
