
#ifndef _libint2_src_bin_libint_vrr1onep1_h_
#define _libint2_src_bin_libint_vrr1onep1_h_

#include <generic_rr.h>
#include <onep_1_1.h>

using namespace std;

namespace libint2 {

  /** VRR Recurrence Relation for 1-e overlap integrals.
    * \tparam where specifies whether angular momentum is decreased, in bra or ket.
  */
    template <class BFSet, FunctionPosition where>
      class VRR_1_Overlap_1 : public GenericRecurrenceRelation< VRR_1_Overlap_1<BFSet,where>,
                                                                BFSet,
                                                                GenIntegralSet_1_1<BFSet,OverlapOper,EmptySet> >
    {
    public:
      typedef VRR_1_Overlap_1 ThisType;
      typedef BFSet BasisFunctionType;
      typedef GenIntegralSet_1_1<BFSet,OverlapOper,EmptySet> TargetType;
      typedef GenericRecurrenceRelation<ThisType,BFSet,TargetType> ParentType;
      friend class GenericRecurrenceRelation<ThisType,BFSet,TargetType>;
      static const unsigned int max_nchildren = 9;

      using ParentType::Instance;

      /// Default directionality
      static bool directional() { return ParentType::default_directional(); }

    private:
      using ParentType::RecurrenceRelation::expr_;
      using ParentType::RecurrenceRelation::nflops_;
      using ParentType::target_;
      using ParentType::is_simple;

      /// Constructor is private, used by ParentType::Instance that mainains registry of these objects
      VRR_1_Overlap_1(const SafePtr<TargetType>&, unsigned int dir);

      static std::string descr() { return "OSVRROverlap"; }
      /** Re-Implementation of GenericRecurrenceRelation::generate_label():
          TwoPRep VRR recurrence relations codes are independent of m (it never appears anywhere in equations), hence
          to avoid generating identical code make sure that the (unique) label has m=0. */
      std::string generate_label() const
      {
        static SafePtr<EmptySet> aux0(new EmptySet);
        ostringstream os;
        os << descr() << to_string(where)
           << genintegralset_label(target_->bra(),target_->ket(),aux0,target_->oper());
        return os.str();
      }
    };

  template <class F, FunctionPosition where>
  VRR_1_Overlap_1<F,where>::VRR_1_Overlap_1(const SafePtr< TargetType >& Tint,
                                            unsigned int dir) :
    ParentType(Tint,dir)
    {
      using namespace libint2::algebra;
      using namespace libint2::prefactor;
      using namespace libint2::braket;
      const F& _1 = unit<F>(dir);

      { // can't apply to contracted basis functions
        F a(Tint->bra(0,0));
        F b(Tint->ket(0,0));
        if (a.contracted() ||
            b.contracted())
          return;
      }

      // if derivative integrals, there will be extra terms (Eq. (143) in Obara & Saika JCP 89)
      const OriginDerivative dA = Tint->bra(0,0).deriv();
      const OriginDerivative dB = Tint->ket(0,0).deriv();
      const bool deriv = dA.zero() == false ||
          dB.zero() == false;

      typedef TargetType ChildType;
      ChildFactory<ThisType,ChildType> factory(this);

      // Build on A or B
      {
        // bf quantum on the build center subtracted by 1
        auto a = ( where == InBra ? Tint->bra(0,0) - _1 : Tint->bra(0,0) );
        if (!exists(a)) return;
        auto b = ( where == InKet ? Tint->ket(0,0) - _1 : Tint->ket(0,0) );
        if (!exists(b)) return;

        auto AB = factory.make_child(a,b);
        if (is_simple()) { expr_ = Vector(where == InBra ? "PA" : "PB")[dir] * AB; nflops_+=1; }

        auto am1 = a - _1;
        if (exists(am1)) {
          auto Am1B = factory.make_child(am1,b);
          if (is_simple()) { expr_ += Vector(a)[dir] * Scalar("oo2z") * Am1B;  nflops_+=3; }
        }
        auto bm1 = b - _1;
        if (exists(bm1)) {
          auto ABm1 = factory.make_child(a,bm1);
          if (is_simple()) { expr_ += Vector(b)[dir] * Scalar("oo2z") * ABm1;  nflops_+=3; }
        }
      }

      // if got here, can decrement by at least 1 quantum
      // add additional derivative terms
      if (deriv) {
        // bf quantum on the build center subtracted by 1
        F a( where == InBra ? Tint->bra(0,0) - _1 : Tint->bra(0,0) );
        F b( where == InKet ? Tint->ket(0,0) - _1 : Tint->ket(0,0) );

        // treatment of derivative terms differs for shell sets and integrals
        // since in computing shell sets transfer/build will occur in all 3 directions
        // change in up to all three derivative indices will occur
        for(unsigned int dxyz=0; dxyz<3; ++dxyz) {

          if (is_simple() && dxyz != dir) // for integrals only consider derivatives in THE build direction
            continue;

          OriginDerivative _d1; _d1.inc(dxyz);

          SafePtr<DGVertex> _nullptr;

          // dA - _1?
          {
            const OriginDerivative dAm1(dA - _d1);
            if (exists(dAm1)) { // yes
              a.deriv() = dAm1;
              auto AB = factory.make_child(a,b);
              if (is_simple()) {
                if (where == InBra) { // building on A
                  expr_ -= Vector(dA)[dxyz] * Scalar("rho12_over_alpha1") * AB;  nflops_ += 3; }
                if (where == InKet) { // building on B
                  expr_ += Vector(dA)[dxyz] * Scalar("rho12_over_alpha2") * AB;  nflops_ += 3; }
              }
              a.deriv() = dA;
            }
          }

          // dB - _1?
          {
            const OriginDerivative dBm1(dB - _d1);
            if (exists(dBm1)) { // yes
              b.deriv() = dBm1;
              auto AB = factory.make_child(a,b);
              if (is_simple()) {
                if (where == InBra) { // building on A
                  expr_ += Vector(dB)[dxyz] * Scalar("rho12_over_alpha1") * AB;  nflops_ += 3; }
                if (where == InKet) { // building on B
                  expr_ -= Vector(dB)[dxyz] * Scalar("rho12_over_alpha2") * AB;  nflops_ += 3; }
              }
              b.deriv() = dB;
            }
          }

        }
      } // end of deriv

      return;
    }


  /** VRR Recurrence Relation for 1-e kinetic energy integrals.
    * \tparam where specifies whether angular momentum is decreased, in bra or ket.
  */
    template <class BFSet, FunctionPosition where>
      class VRR_1_Kinetic_1 : public GenericRecurrenceRelation< VRR_1_Kinetic_1<BFSet,where>,
                                                                BFSet,
                                                                GenIntegralSet_1_1<BFSet,KineticOper,EmptySet> >
    {
    public:
      typedef VRR_1_Kinetic_1 ThisType;
      typedef BFSet BasisFunctionType;
      typedef GenIntegralSet_1_1<BFSet,KineticOper,EmptySet> TargetType;
      typedef GenericRecurrenceRelation<ThisType,BFSet,TargetType> ParentType;
      friend class GenericRecurrenceRelation<ThisType,BFSet,TargetType>;
      static const unsigned int max_nchildren = 9;

      using ParentType::Instance;

      /// Default directionality
      static bool directional() { return ParentType::default_directional(); }

    private:
      using ParentType::RecurrenceRelation::expr_;
      using ParentType::RecurrenceRelation::nflops_;
      using ParentType::target_;
      using ParentType::is_simple;

      /// Constructor is private, used by ParentType::Instance that mainains registry of these objects
      VRR_1_Kinetic_1(const SafePtr<TargetType>&, unsigned int dir);

      static std::string descr() { return "OSVRRKinetic"; }
      /** Re-Implementation of GenericRecurrenceRelation::generate_label():
          TwoPRep VRR recurrence relations codes are independent of m (it never appears anywhere in equations), hence
          to avoid generating identical code make sure that the (unique) label has m=0. */
      std::string generate_label() const
      {
        static SafePtr<EmptySet> aux0(new EmptySet);
        ostringstream os;
        os << descr() << to_string(where)
           << genintegralset_label(target_->bra(),target_->ket(),aux0,target_->oper());
        return os.str();
      }
    };

  template <class F, FunctionPosition where>
  VRR_1_Kinetic_1<F,where>::VRR_1_Kinetic_1(const SafePtr< TargetType >& Tint,
                                            unsigned int dir) :
    ParentType(Tint,dir)
    {
      using namespace libint2::algebra;
      using namespace libint2::prefactor;
      using namespace libint2::braket;
      const F& _1 = unit<F>(dir);

      { // can't apply to contracted basis functions
        F a(Tint->bra(0,0));
        F b(Tint->ket(0,0));
        if (a.contracted() ||
            b.contracted())
          return;
      }

      // if derivative integrals, there will be extra terms (Eq. (143) in Obara & Saika JCP 89)
      const OriginDerivative dA = Tint->bra(0,0).deriv();
      const OriginDerivative dB = Tint->ket(0,0).deriv();
      const bool deriv = dA.zero() == false ||
          dB.zero() == false;

      typedef TargetType Child1Type;
      ChildFactory<ThisType,Child1Type> factory(this);
      typedef GenIntegralSet_1_1<F,OverlapOper,EmptySet> Child2Type;
      ChildFactory<ThisType,Child2Type> overlap_factory(this);

      // Build on A or B
      {
        // bf quantum on the build center subtracted by 1
        auto a = ( where == InBra ? Tint->bra(0,0) - _1 : Tint->bra(0,0) );
        if (!exists(a)) return;
        auto b = ( where == InKet ? Tint->ket(0,0) - _1 : Tint->ket(0,0) );
        if (!exists(b)) return;

        auto AB = factory.make_child(a,b);
        if (is_simple()) { expr_ = Vector(where == InBra ? "PA" : "PB")[dir] * AB; nflops_+=1; }

        auto am1 = a - _1;
        if (exists(am1)) {
          auto Am1B = factory.make_child(am1,b);
          auto S_Am1B = (where == InBra) ? overlap_factory.make_child(am1,b) : 0;
          if (is_simple()) {
            if (where == InBra) {
              expr_ += Vector(a)[dir] * ( Scalar("oo2z") * Am1B -
                                          Scalar("rho12_over_alpha1") * S_Am1B );
              nflops_+=5;
            }
            else {
              expr_ += Vector(a)[dir] * Scalar("oo2z") * Am1B;
              nflops_+=3;
            }
          }
        }
        auto bm1 = b - _1;
        if (exists(bm1)) {
          auto ABm1 = factory.make_child(a,bm1);
          auto S_ABm1 = (where == InKet) ? overlap_factory.make_child(a,bm1) : 0;
          if (is_simple()) {
            if (where == InKet) {
              expr_ += Vector(b)[dir] * ( Scalar("oo2z") * ABm1 -
                                          Scalar("rho12_over_alpha2") * S_ABm1 );
              nflops_+=5;
            }
            else {
              expr_ += Vector(b)[dir] * Scalar("oo2z") * ABm1;
              nflops_+=3;
            }
          }
        }

        {
          auto S_AB_target = where == InBra ? overlap_factory.make_child(a + _1,b) : overlap_factory.make_child(a,b + _1);
          if (is_simple()) {
            expr_ += Scalar("two_rho12") * S_AB_target;
            nflops_+=2;
          }
        }

      }

      // if got here, can decrement by at least 1 quantum
      // add additional derivative terms
      if (deriv) {
        assert(false); // not yet implemented

        // bf quantum on the build center subtracted by 1
        F a( where == InBra ? Tint->bra(0,0) - _1 : Tint->bra(0,0) );
        F b( where == InKet ? Tint->ket(0,0) - _1 : Tint->ket(0,0) );

        // treatment of derivative terms differs for shell sets and integrals
        // since in computing shell sets transfer/build will occur in all 3 directions
        // change in up to all three derivative indices will occur
        for(unsigned int dxyz=0; dxyz<3; ++dxyz) {

          if (is_simple() && dxyz != dir) // for integrals only consider derivatives in THE build direction
            continue;

          OriginDerivative _d1; _d1.inc(dxyz);

          SafePtr<DGVertex> _nullptr;

          // dA - _1?
          {
            const OriginDerivative dAm1(dA - _d1);
            if (exists(dAm1)) { // yes
              a.deriv() = dAm1;
              auto AB = factory.make_child(a,b);
              if (is_simple()) {
                if (where == InBra) { // building on A
                  expr_ -= Vector(dA)[dxyz] * Scalar("rho12_over_alpha1") * AB;  nflops_ += 3; }
                if (where == InKet) { // building on B
                  expr_ += Vector(dA)[dxyz] * Scalar("rho12_over_alpha2") * AB;  nflops_ += 3; }
              }
              a.deriv() = dA;
            }
          }

          // dB - _1?
          {
            const OriginDerivative dBm1(dB - _d1);
            if (exists(dBm1)) { // yes
              b.deriv() = dBm1;
              auto AB = factory.make_child(a,b);
              if (is_simple()) {
                if (where == InBra) { // building on A
                  expr_ += Vector(dB)[dxyz] * Scalar("rho12_over_alpha1") * AB;  nflops_ += 3; }
                if (where == InKet) { // building on B
                  expr_ -= Vector(dB)[dxyz] * Scalar("rho12_over_alpha2") * AB;  nflops_ += 3; }
              }
              b.deriv() = dB;
            }
          }

        }
      } // end of deriv

      return;
    }

}; // namespace libint2

#endif
