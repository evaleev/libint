/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

#ifndef _libint2_src_bin_libint_prefactors_h_
#define _libint2_src_bin_libint_prefactors_h_

#include <cstring>
#include <smart_ptr.h>
#include <entity.h>
#include <bfset.h>
#include <singl_stack.h>

#define CTIMEENTITIES_SINGLETONS 1

namespace libint2 {

  /**
     Prefactors is a collection of common quantities which appear
     as prefactors in recurrence relations for Gaussian integrals.
     See Obara-Saika paper for description of the most common ones.
  */

  class Prefactors {
  public:
    Prefactors();
    ~Prefactors();

    typedef RTimeEntity<double> rdouble;
    typedef CTimeEntity<double> cdouble;

    static const unsigned int NMAX = 200;
    static const unsigned int np = 2;
    static const unsigned int nfunc_per_part = 1;

    /**
       X-Y vectors, where X and Y are for the same particle:
       X_Y[0] = AB, X_Y[1] = CD, etc.
     */
    SafePtr<rdouble> vX_Y[np];
    /// Cartesian components of X-Y vectors
    SafePtr<rdouble> X_Y[np][3];

    /**
       Y-X vectors, where X and Y are for the same particle:
       Y_X[0] = BA, Y_X[1] = DC, etc.
     */
    SafePtr<rdouble> vY_X[np];
    /// Cartesian components of Y_X vectors
    SafePtr<rdouble> Y_X[np][3];

    /**
       XY-X vectors:
       XY is either P or Q,
       X is either (A or B) or (C or D).
       Hence, vXY_X[0][0] = P-A, vXY_X[1][0] = Q-C,
       vXY_X[0][1] = P-B, and vXY_X[1][1] = Q-D.
    */
    SafePtr<rdouble> vXY_X[np][2];
    /// cartesian components of vXY_X vector
    SafePtr<rdouble> XY_X[np][2][3];
    /**
       W-XY vectors:
       vW_XY[0] = W-P, vW_XY[1] = W-Q.
    */
    SafePtr<rdouble> vW_XY[np];
    /// cartesian components of W_XY vector
    SafePtr<rdouble> W_XY[np][3];

    /**
       orbital exponents
    */
    SafePtr<rdouble> zeta[np][2];
    /**
       squared orbital exponents
    */
    SafePtr<rdouble> zeta2[np][2];

    /**
       alpha12[p] is the sum of exponents for particle p:
       alpha12[0] = zeta,
       alpha12[1] = eta.
    */
    SafePtr<rdouble> alpha12[np];
    /// rho = zeta*eta/(zeta+eta)
    SafePtr<rdouble> rho;
    /// 1/(2*alpha12)
    SafePtr<rdouble> one_o_2alpha12[np];
    /// rho/alpha12
    SafePtr<rdouble> rho_o_alpha12[np];
    /// 1/(2*(zeta+eta))
    SafePtr<rdouble> one_o_2alphasum;

    /**
    Prefactors for the ITR relation for TwoPRep integrals (a+1 0|c0):
    */
    /// prefactor in front of (a0|c0) = -(zeta[0][1] AB + zeta[1][1] CD)/alpha12[0]
    SafePtr<rdouble> TwoPRepITR_vpfac0[np];
    /// cartesian components of pfac0 vector
    SafePtr<rdouble> TwoPRepITR_pfac0[np][3];
    /// prefactor in front of (a0|c+1 0) = -alpha12[1]/alpha12[0]
    SafePtr<rdouble> TwoPRepITR_pfac1[np];
    /// prefactor in front of (a-1 0|c0) is one_o_2alpha12[0]
    /// prefactor in front of (a0|c-1 0) is one_o_2alpha12[0]

    /**
       Prefactors for the VRR relation for R12_k_G12 integrals (k>=0):
    */
    /// prefactor in front of (a0|c0)
    SafePtr<rdouble> R12kG12VRR_vpfac0[np];
    /// cartesian components of pfac0 vector
    SafePtr<rdouble> R12kG12VRR_pfac0[np][3];
    /// prefactor in front of (a-1 0|c0)
    SafePtr<rdouble> R12kG12VRR_pfac1[np];
    /// prefactor in front of (a0|c-1 0)
    SafePtr<rdouble> R12kG12VRR_pfac2;
    /// prefactor in front of (|k-2|)
    SafePtr<rdouble> R12kG12VRR_pfac3[np];
    /// prefactor in front of (a0|k-2|c0)
    SafePtr<rdouble> R12kG12VRR_vpfac4[np];
    /// cartesian components of pfac4 vector
    SafePtr<rdouble> R12kG12VRR_pfac4[np][3];

    /**
     * Precomputed 1-d integrals
     */
    /// (0|0)_xyz 1-d overlap integrals
    SafePtr<rdouble> Overlap00_1d[3];


#if CTIMEENTITIES_SINGLETONS
    /// integers represented as doubles
    SafePtr<cdouble> N_i[NMAX];

    SafePtr<cdouble> Cdouble(double a);
#endif

  private:

  };

  extern Prefactors prefactors;

  namespace prefactor {

    template <typename T> struct RTimeSingletons {
      typedef SingletonStack<RTimeEntity<T>,typename RTimeEntity<T>::key_type> ManagerType;
      static SafePtr<ManagerType>& Manager() {
        if (manager_ == 0) {
          manager_ = SafePtr<ManagerType>(new ManagerType(&RTimeEntity<T>::key));
        }
        return manager_;
      }
      static SafePtr<ManagerType> manager_;
    };
    template <typename T> SafePtr<typename RTimeSingletons<T>::ManagerType> RTimeSingletons<T>::manager_;

#if CTIMEENTITIES_SINGLETONS
    template <typename T> struct CTimeSingletons {
      typedef SingletonStack<CTimeEntity<T>,T> ManagerType;
      static SafePtr<ManagerType>& Manager() {
        if (manager_ == 0) {
          manager_ = SafePtr<ManagerType>(new ManagerType(&CTimeEntity<T>::value));
        }
        return manager_;
      }
      static SafePtr<ManagerType> manager_;
    };
    template <typename T> SafePtr<typename CTimeSingletons<T>::ManagerType> CTimeSingletons<T>::manager_;
#endif

    /// make a floating-point compile-time quantity from an integer
    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value>::type>
    SafePtr< CTimeEntity<double> > Scalar(T a) {
      typedef CTimeEntity<double> return_type;
      SafePtr<return_type> tmp(new return_type(a));
#if CTIMEENTITIES_SINGLETONS
      typedef CTimeSingletons<double> singletons_type;
      typedef typename singletons_type::ManagerType ManagerType;
      const typename ManagerType::value_type& result = singletons_type::Manager()->find(tmp);
      return result.second;
#else
      return tmp;
#endif
    }
    /// make a compile-time quantity
    inline SafePtr< CTimeEntity<double> > Scalar(double a) {
      typedef CTimeEntity<double> return_type;
      SafePtr<return_type> tmp(new return_type(a));
#if CTIMEENTITIES_SINGLETONS
      typedef CTimeSingletons<double> singletons_type;
      typedef typename singletons_type::ManagerType ManagerType;
      const typename ManagerType::value_type& result = singletons_type::Manager()->find(tmp);
      return result.second;
#else
      return tmp;
#endif
    }
    /// make a runtime quantity
//    inline SafePtr< RTimeEntity<double> > Scalar(const char* id) {
//      typedef double T;
//      typedef RTimeEntity<T> return_type;
//      typedef RTimeSingletons<T> singletons_type;
//      typedef singletons_type::ManagerType ManagerType;
//      SafePtr<return_type> tmp(new return_type(id));
//      const ManagerType::value_type& result = singletons_type::Manager()->find(tmp);
//      return result.second;
//    }
    /// make a runtime quantity
    inline SafePtr< RTimeEntity<double> > Scalar(const std::string& id) {
      typedef double T;
      typedef RTimeEntity<T> return_type;
      typedef RTimeSingletons<T> singletons_type;
      typedef singletons_type::ManagerType ManagerType;
      SafePtr<return_type> tmp(new return_type(id));
      const ManagerType::value_type& result = singletons_type::Manager()->find(tmp);
      return result.second;
    }

    /// auxiliary class that write expressions with runtime cartesian vectors
    template <class T> class RTimeVector3 {
    public:
      RTimeVector3(const char* id) : id_(id) { }
      RTimeVector3(const std::string& id) : id_(id) { }
      SafePtr< RTimeEntity<T> > operator[](unsigned int xyz) {
        return Scalar(id_ + "_" + dirchar[xyz]);
      }
      private:
        static const char* dirchar;
        std::string id_;
    };
    template <class T> const char* RTimeVector3<T>::dirchar(strdup("xyz"));

    /// auxiliary class that write expressions with compile-time cartesian vectors
    template <class T> class CTimeVector3 {
    public:
      CTimeVector3(const T* val) {
        for(int xyz=0; xyz<3; ++xyz) val_[xyz] = val[xyz];
      }
      SafePtr< CTimeEntity<T> > operator[](unsigned int xyz) {
        return Scalar(val_[xyz]);
      }
      private:
        T val_[3];
    };

    /// make a runtime quantity
#if 0 // these do not get resolved correctly by icpc 12 on OS X
    template <class T = double> RTimeVector3<T> Vector(const char* id)
    {
      return RTimeVector3<T>(id);
    }
    /// make a runtime quantity
    template <class T = double> RTimeVector3<T> Vector(const std::string& id)
    {
      return RTimeVector3<T>(id);
    }
    /// make a compile-time quantity
    template <class T = double> CTimeVector3<T> Vector(const T* a)
    {
      return CTimeVector3<T>(a);
    }
#endif
    inline RTimeVector3<double> Vector(const char* id)
    {
      return RTimeVector3<double>(id);
    }
    /// make a compile-time quantity
    inline CTimeVector3<double> Vector(const CGF& bf)
    {
      double qn[3]; for(unsigned int xyz=0; xyz<3; ++xyz) qn[xyz] = bf.qn(xyz);
      return CTimeVector3<double>(qn);
    }
    /// make a compile-time quantity
    inline CTimeVector3<double> Vector(const OriginDerivative<3u>& dd)
    {
      double d[3]; for(unsigned int xyz=0; xyz<3; ++xyz) d[xyz] = dd.d(xyz);
      return CTimeVector3<double>(d);
    }

  }
};

#endif

