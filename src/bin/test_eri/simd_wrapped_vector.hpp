//Name: Jonathan Dullea 
//jdullea@umass.edu
//This header was generatef by generate_vector.cpp

#include <iostream>
#include <x86intrin.h>
#include "immintrin.h"
#include <cstring>
#include <cmath>
using namespace std;

template <typename Real, unsigned int n>
class VectorSIMD{
	VectorSIMD(){};
	VectorSIMD(Real a) {}};

template <typename Real>
class VectorSIMD<Real,2>{
public:
		__m256d _avx0;
		VectorSIMD(){};
		VectorSIMD(Real a){
			_avx0=_mm256_set_pd(0,0,a,a);
		}

		VectorSIMD(Real (&a)[2]){
			_avx0=_mm256_set_pd(0,0,a[1],a[0]);
		}

		VectorSIMD(Real _0,Real _1){
			_avx0=_mm256_set_pd( 0,0,_1,_0);
		}

		VectorSIMD(__m256d _0){
			_avx0= _0;
		}

		VectorSIMD& operator=(Real a){
			_avx0 = _mm256_set_pd(0,0,a,a);
			return *this;
		}

		VectorSIMD& operator += (VectorSIMD<Real,2> a){
			_avx0=  _mm256_add_pd(_avx0,a._avx0);
			return *this;
		}

		VectorSIMD& operator -= (VectorSIMD<Real,2> a){
			_avx0=  _mm256_sub_pd(_avx0,a._avx0);
			return *this;
		}

		VectorSIMD operator -() const{
			const static __m256d m1 = _mm256_set_pd(-1.0,-1.0,-1.0,-1.0);
			VectorSIMD<double,2> result;
			result._avx0=  _mm256_mul_pd(this->_avx0,m1);
			return result;
		}

		void load(Real const* a){
			_avx0 =  _mm256_loadu_pd(&a[0]);
		}

		void load_aligned(Real const* a){
			_avx0 =  _mm256_load_pd(&a[0]);
		}

		void convert(Real *a) const {
			double temp[4];
			_mm256_storeu_pd(&temp[0],_avx0);
			a[0] = temp[0];
			a[1] = temp[1];
		}

		void convert_aligned(double *a) const {
			double temp[4];
			_mm256_storeu_pd(&temp[0],_avx0);
			a[0] = temp[0];
			a[1] = temp[1];
		}

	};


	inline VectorSIMD<double,2> operator*(double a, VectorSIMD<double,2> b){
		VectorSIMD<double,2> c;
		__m256d _a =  _mm256_set_pd(a,a,a,a);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		return c;
	}


	inline VectorSIMD<double,2> operator*(VectorSIMD<double,2> a ,double b){
		VectorSIMD<double,2> c;
		__m256d _b =  _mm256_set_pd(b,b,b,b);
		c._avx0=  _mm256_mul_pd(a._avx0, _b);
		return c;
	}


	inline VectorSIMD<double,2> operator*(int a, VectorSIMD<double,2> b){
		if(a==1){return b;}
		VectorSIMD<double,2> c;
		double q = static_cast<double>(a);
		__m256d _a =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		return c;
	}


	inline VectorSIMD<double,2> operator*(VectorSIMD<double,2> a, int b){
		if(b==1){return a;}
		VectorSIMD<double,2> c;
		double q = static_cast<double>(b);
		__m256d _b =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_b, a._avx0);
		return c;
	}


	inline VectorSIMD<double,2> operator*(VectorSIMD<double,2> a, VectorSIMD<double,2> b){
		VectorSIMD<double,2> c;
		c._avx0=  _mm256_mul_pd(a._avx0, b._avx0);
		return c;
	}


	inline VectorSIMD<double,2> operator+(VectorSIMD<double,2> a, VectorSIMD<double,2> b){
		VectorSIMD<double,2> c;
		c._avx0=  _mm256_add_pd(a._avx0, b._avx0);
		return c;
	}


	inline VectorSIMD<double,2> operator-(VectorSIMD<double,2> a, VectorSIMD<double,2> b){
		VectorSIMD<double,2> c;
		c._avx0=  _mm256_sub_pd(a._avx0, b._avx0);
		return c;
	}


	inline VectorSIMD<double,2> operator/(VectorSIMD<double,2> a, VectorSIMD<double,2> b){
		VectorSIMD<double,2> c;
		c._avx0=  _mm256_div_pd(a._avx0, b._avx0);
		return c;
	}


	inline VectorSIMD<double,2> exp(VectorSIMD<double,2> a){
		double a_d[2]; a.convert(a_d);
		a_d[0] = std::exp(a_d[0]);
		a_d[1] = std::exp(a_d[1]);
		VectorSIMD<double,2> r(a_d);
		return r;
	}

	inline VectorSIMD<double,2> sqrt(VectorSIMD<double,2> a){
		double a_d[2]; a.convert(a_d);
		a_d[0] = std::sqrt(a_d[0]);
		a_d[1] = std::sqrt(a_d[1]);
		VectorSIMD<double,2> r(a_d);
		return r;
	}

	inline std::ostream& operator<<(std::ostream& os, VectorSIMD<double,2> a){ 
		double ad[2];
		a.convert(ad);
		os << "{" << ad[0] << ","  << ad[1] << "}";
		return os;
		}
	
	inline double reduce(VectorSIMD<double,2> a){
		
	       __m256d s = _mm256_hadd_pd(a._avx0,a._avx0);
	      double new_eri_double =  ((double*)&s)[0] + ((double*)&s)[2];
	      return new_eri_double;
	}	



template <typename Real>
class VectorSIMD<Real,3>{
public:
		__m256d _avx0;
		VectorSIMD(){};
		VectorSIMD(Real a){
			_avx0=_mm256_set_pd(0,a,a,a);
		}

		VectorSIMD(Real (&a)[3]){
			_avx0=_mm256_set_pd(0,a[2],a[1],a[0]);
		}

		VectorSIMD(Real _0,Real _1,Real _2){
			_avx0=_mm256_set_pd( 0,_2,_1,_0);
		}

		VectorSIMD(__m256d _0){
			_avx0= _0;
		}

		VectorSIMD& operator=(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			return *this;
		}

		VectorSIMD& operator += (VectorSIMD<Real,3> a){
			_avx0=  _mm256_add_pd(_avx0,a._avx0);
			return *this;
		}

		VectorSIMD& operator -= (VectorSIMD<Real,3> a){
			_avx0=  _mm256_sub_pd(_avx0,a._avx0);
			return *this;
		}

		VectorSIMD operator -() const{
			const static __m256d m1 = _mm256_set_pd(-1.0,-1.0,-1.0,-1.0);
			VectorSIMD<double,3> result;
			result._avx0=  _mm256_mul_pd(this->_avx0,m1);
			return result;
		}

		void load(Real const* a){
			_avx0 =  _mm256_loadu_pd(&a[0]);
		}

		void load_aligned(Real const* a){
			_avx0 =  _mm256_load_pd(&a[0]);
		}

		void convert(Real *a) const {
			_mm256_storeu_pd(&a[0],_avx0);
		}

		void convert_aligned(Real *a) const {
			_mm256_storeu_pd(&a[0],_avx0);
		}

	};


	inline VectorSIMD<double,3> operator*(double a, VectorSIMD<double,3> b){
		VectorSIMD<double,3> c;
		__m256d _a =  _mm256_set_pd(a,a,a,a);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		return c;
	}


	inline VectorSIMD<double,3> operator*(VectorSIMD<double,3> a ,double b){
		VectorSIMD<double,3> c;
		__m256d _b =  _mm256_set_pd(b,b,b,b);
		c._avx0=  _mm256_mul_pd(a._avx0, _b);
		return c;
	}


	inline VectorSIMD<double,3> operator*(int a, VectorSIMD<double,3> b){
		if(a==1){return b;}
		VectorSIMD<double,3> c;
		double q = static_cast<double>(a);
		__m256d _a =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		return c;
	}


	inline VectorSIMD<double,3> operator*(VectorSIMD<double,3> a, int b){
		if(b==1){return a;}
		VectorSIMD<double,3> c;
		double q = static_cast<double>(b);
		__m256d _b =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_b, a._avx0);
		return c;
	}


	inline VectorSIMD<double,3> operator*(VectorSIMD<double,3> a, VectorSIMD<double,3> b){
		VectorSIMD<double,3> c;
		c._avx0=  _mm256_mul_pd(a._avx0, b._avx0);
		return c;
	}


	inline VectorSIMD<double,3> operator+(VectorSIMD<double,3> a, VectorSIMD<double,3> b){
		VectorSIMD<double,3> c;
		c._avx0=  _mm256_add_pd(a._avx0, b._avx0);
		return c;
	}


	inline VectorSIMD<double,3> operator-(VectorSIMD<double,3> a, VectorSIMD<double,3> b){
		VectorSIMD<double,3> c;
		c._avx0=  _mm256_sub_pd(a._avx0, b._avx0);
		return c;
	}


	inline VectorSIMD<double,3> operator/(VectorSIMD<double,3> a, VectorSIMD<double,3> b){
		VectorSIMD<double,3> c;
		c._avx0=  _mm256_div_pd(a._avx0, b._avx0);
		return c;
	}


	inline VectorSIMD<double,3> exp(VectorSIMD<double,3> a){
		double a_d[3]; a.convert(a_d);
		a_d[0] = std::exp(a_d[0]);
		a_d[1] = std::exp(a_d[1]);
		a_d[2] = std::exp(a_d[2]);
		VectorSIMD<double,3> r(a_d);
		return r;
	}

	inline VectorSIMD<double,3> sqrt(VectorSIMD<double,3> a){
		double a_d[3]; a.convert(a_d);
		a_d[0] = std::sqrt(a_d[0]);
		a_d[1] = std::sqrt(a_d[1]);
		a_d[2] = std::sqrt(a_d[2]);
		VectorSIMD<double,3> r(a_d);
		return r;
	}

	inline std::ostream& operator<<(std::ostream& os, VectorSIMD<double,3> a){ 
		double ad[3];
		a.convert(ad);
		os << "{" << ad[0] << ","  << ad[1] << ","  << ad[2] << "}";
		return os;
		}


	inline double reduce(VectorSIMD<double,3> a){
	       __m256d s = _mm256_hadd_pd(a._avx0,a._avx0);
	      double new_eri_double =  ((double*)&s)[0] + ((double*)&s)[2];
	      return new_eri_double;
	}	



template <typename Real>
class VectorSIMD<Real,4>{
public:
		__m256d _avx0;
		VectorSIMD(){};
		VectorSIMD(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
		}

		VectorSIMD(Real (&a)[4]){
			_avx0=_mm256_set_pd(a[3],a[2],a[1],a[0]);
		}

		VectorSIMD(Real _0,Real _1,Real _2,Real _3){
			_avx0=_mm256_set_pd(_3,_2,_1,_0);
		}

		VectorSIMD(__m256d _0){
			_avx0= _0;
		}

		VectorSIMD& operator=(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			return *this;
		}

		VectorSIMD& operator += (VectorSIMD<Real,4> a){
			_avx0=  _mm256_add_pd(_avx0,a._avx0);
			return *this;
		}

		VectorSIMD& operator -= (VectorSIMD<Real,4> a){
			_avx0=  _mm256_sub_pd(_avx0,a._avx0);
			return *this;
		}

		VectorSIMD operator -() const{
			const static __m256d m1 = _mm256_set_pd(-1.0,-1.0,-1.0,-1.0);
			VectorSIMD<double,4> result;
			result._avx0=  _mm256_mul_pd(this->_avx0,m1);
			return result;
		}

		void load(Real const* a){
			_avx0 =  _mm256_loadu_pd(&a[0]);
		}

		void load_aligned(Real const* a){
			_avx0 =  _mm256_load_pd(&a[0]);
		}

		void convert(Real *a) const {
			_mm256_storeu_pd(&a[0],_avx0);
		}

		void convert_aligned(Real *a) const {
			_mm256_storeu_pd(&a[0],_avx0);
		}

	};


	inline VectorSIMD<double,4> operator*(double a, VectorSIMD<double,4> b){
		VectorSIMD<double,4> c;
		__m256d _a =  _mm256_set_pd(a,a,a,a);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		return c;
	}


	inline VectorSIMD<double,4> operator*(VectorSIMD<double,4> a ,double b){
		VectorSIMD<double,4> c;
		__m256d _b =  _mm256_set_pd(b,b,b,b);
		c._avx0=  _mm256_mul_pd(a._avx0, _b);
		return c;
	}


	inline VectorSIMD<double,4> operator*(int a, VectorSIMD<double,4> b){
		if(a==1){return b;}
		VectorSIMD<double,4> c;
		double q = static_cast<double>(a);
		__m256d _a =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		return c;
	}


	inline VectorSIMD<double,4> operator*(VectorSIMD<double,4> a, int b){
		if(b==1){return a;}
		VectorSIMD<double,4> c;
		double q = static_cast<double>(b);
		__m256d _b =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_b, a._avx0);
		return c;
	}


	inline VectorSIMD<double,4> operator*(VectorSIMD<double,4> a, VectorSIMD<double,4> b){
		VectorSIMD<double,4> c;
		c._avx0=  _mm256_mul_pd(a._avx0, b._avx0);
		return c;
	}


	inline VectorSIMD<double,4> operator+(VectorSIMD<double,4> a, VectorSIMD<double,4> b){
		VectorSIMD<double,4> c;
		c._avx0=  _mm256_add_pd(a._avx0, b._avx0);
		return c;
	}


	inline VectorSIMD<double,4> operator-(VectorSIMD<double,4> a, VectorSIMD<double,4> b){
		VectorSIMD<double,4> c;
		c._avx0=  _mm256_sub_pd(a._avx0, b._avx0);
		return c;
	}


	inline VectorSIMD<double,4> operator/(VectorSIMD<double,4> a, VectorSIMD<double,4> b){
		VectorSIMD<double,4> c;
		c._avx0=  _mm256_div_pd(a._avx0, b._avx0);
		return c;
	}


	inline VectorSIMD<double,4> exp(VectorSIMD<double,4> a){
		double a_d[4]; a.convert(a_d);
		a_d[0] = std::exp(a_d[0]);
		a_d[1] = std::exp(a_d[1]);
		a_d[2] = std::exp(a_d[2]);
		a_d[3] = std::exp(a_d[3]);
		VectorSIMD<double,4> r(a_d);
		return r;
	}

	inline VectorSIMD<double,4> sqrt(VectorSIMD<double,4> a){
		double a_d[4]; a.convert(a_d);
		a_d[0] = std::sqrt(a_d[0]);
		a_d[1] = std::sqrt(a_d[1]);
		a_d[2] = std::sqrt(a_d[2]);
		a_d[3] = std::sqrt(a_d[3]);
		VectorSIMD<double,4> r(a_d);
		return r;
	}

	inline std::ostream& operator<<(std::ostream& os, VectorSIMD<double,4> a){ 
		double ad[4];
		a.convert(ad);
		os << "{" << ad[0] << ","  << ad[1] << ","  << ad[2] << ","  << ad[3] << "}";
		return os;
	}

	inline double reduce(VectorSIMD<double,4> a){
	       __m256d s = _mm256_hadd_pd(a._avx0,a._avx0);
	      double new_eri_double =  ((double*)&s)[0] + ((double*)&s)[2];
	      return new_eri_double;
	}	



template <typename Real>
class VectorSIMD<Real,5>{
public:
		__m256d _avx0;
		__m256d _avx1;
		VectorSIMD(){};
		VectorSIMD(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			_avx1=_mm256_set_pd(0,0,0,a);
		}

		VectorSIMD(Real (&a)[5]){
			_avx0=_mm256_set_pd(a[3],a[2],a[1],a[0]);
			_avx1=_mm256_set_pd(0,0,0,a[4]);
		}

		VectorSIMD(Real _0,Real _1,Real _2,Real _3,Real _4){
			_avx0=_mm256_set_pd(_3,_2,_1,_0);
			_avx1=_mm256_set_pd( 0,0,0,_4);
		}

		VectorSIMD(__m256d _0,__m256d _1){
			_avx0= _0;
			_avx1= _1;
		}

		VectorSIMD& operator=(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			_avx1 = _mm256_set_pd(0,0,0,a);
			return *this;
		}

		VectorSIMD& operator += (VectorSIMD<Real,5> a){
			_avx0=  _mm256_add_pd(_avx0,a._avx0);
			_avx1=  _mm256_add_pd(_avx1,a._avx1);
			return *this;
		}

		VectorSIMD& operator -= (VectorSIMD<Real,5> a){
			_avx0=  _mm256_sub_pd(_avx0,a._avx0);
			_avx1=  _mm256_sub_pd(_avx1,a._avx1);
			return *this;
		}

		VectorSIMD operator -() const{
			const static __m256d m1 = _mm256_set_pd(-1.0,-1.0,-1.0,-1.0);
			VectorSIMD<double,5> result;
			result._avx0=  _mm256_mul_pd(this->_avx0,m1);
			result._avx1=  _mm256_mul_pd(this->_avx1,m1);
			return result;
		}

		void load(Real const* a){
			_avx0 =  _mm256_loadu_pd(&a[0]);
			_avx1 =  _mm256_loadu_pd(&a[4]);
		}

		void load_aligned(Real const* a){
			_avx0 =  _mm256_load_pd(&a[0]);
			_avx1 =  _mm256_load_pd(&a[4]);
		}

		void convert(Real *a) const {
			double temp[8];
			_mm256_storeu_pd(&temp[0],_avx0);
			_mm256_storeu_pd(&temp[4],_avx1);
			a[0] = temp[0];
			a[1] = temp[1];
			a[2] = temp[2];
			a[3] = temp[3];
			a[4] = temp[4];
		}

		void convert_aligned(double *a) const {
			double temp[8];
			_mm256_storeu_pd(&temp[0],_avx0);
			_mm256_storeu_pd(&temp[4],_avx1);
			a[0] = temp[0];
			a[1] = temp[1];
			a[2] = temp[2];
			a[3] = temp[3];
			a[4] = temp[4];
		}

	};


	inline VectorSIMD<double,5> operator*(double a, VectorSIMD<double,5> b){
		VectorSIMD<double,5> c;
		__m256d _a =  _mm256_set_pd(a,a,a,a);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		c._avx1=  _mm256_mul_pd(_a, b._avx1);
		return c;
	}


	inline VectorSIMD<double,5> operator*(VectorSIMD<double,5> a ,double b){
		VectorSIMD<double,5> c;
		__m256d _b =  _mm256_set_pd(b,b,b,b);
		c._avx0=  _mm256_mul_pd(a._avx0, _b);
		c._avx1=  _mm256_mul_pd(a._avx1, _b);
		return c;
	}


	inline VectorSIMD<double,5> operator*(int a, VectorSIMD<double,5> b){
		if(a==1){return b;}
		VectorSIMD<double,5> c;
		double q = static_cast<double>(a);
		__m256d _a =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		c._avx1=  _mm256_mul_pd(_a, b._avx1);
		return c;
	}


	inline VectorSIMD<double,5> operator*(VectorSIMD<double,5> a, int b){
		if(b==1){return a;}
		VectorSIMD<double,5> c;
		double q = static_cast<double>(b);
		__m256d _b =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_b, a._avx0);
		c._avx1=  _mm256_mul_pd(_b, a._avx1);
		return c;
	}


	inline VectorSIMD<double,5> operator*(VectorSIMD<double,5> a, VectorSIMD<double,5> b){
		VectorSIMD<double,5> c;
		c._avx0=  _mm256_mul_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_mul_pd(a._avx1, b._avx1);
		return c;
	}


	inline VectorSIMD<double,5> operator+(VectorSIMD<double,5> a, VectorSIMD<double,5> b){
		VectorSIMD<double,5> c;
		c._avx0=  _mm256_add_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_add_pd(a._avx1, b._avx1);
		return c;
	}


	inline VectorSIMD<double,5> operator-(VectorSIMD<double,5> a, VectorSIMD<double,5> b){
		VectorSIMD<double,5> c;
		c._avx0=  _mm256_sub_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_sub_pd(a._avx1, b._avx1);
		return c;
	}


	inline VectorSIMD<double,5> operator/(VectorSIMD<double,5> a, VectorSIMD<double,5> b){
		VectorSIMD<double,5> c;
		c._avx0=  _mm256_div_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_div_pd(a._avx1, b._avx1);
		return c;
	}


	inline VectorSIMD<double,5> exp(VectorSIMD<double,5> a){
		double a_d[5]; a.convert(a_d);
		a_d[0] = std::exp(a_d[0]);
		a_d[1] = std::exp(a_d[1]);
		a_d[2] = std::exp(a_d[2]);
		a_d[3] = std::exp(a_d[3]);
		a_d[4] = std::exp(a_d[4]);
		VectorSIMD<double,5> r(a_d);
		return r;
	}

	inline VectorSIMD<double,5> sqrt(VectorSIMD<double,5> a){
		double a_d[5]; a.convert(a_d);
		a_d[0] = std::sqrt(a_d[0]);
		a_d[1] = std::sqrt(a_d[1]);
		a_d[2] = std::sqrt(a_d[2]);
		a_d[3] = std::sqrt(a_d[3]);
		a_d[4] = std::sqrt(a_d[4]);
		VectorSIMD<double,5> r(a_d);
		return r;
	}

	inline std::ostream& operator<<(std::ostream& os, VectorSIMD<double,5> a){ 
		double ad[5];
		a.convert(ad);
		os << "{" << ad[0] << ","  << ad[1] << ","  << ad[2] << ","  << ad[3] << ","  << ad[4] << "}";
		return os;
		}

	inline double reduce(VectorSIMD<double,5> a){
		__m256d mask = _mm256_set_pd(0,0,0,1); //this ensures that assumed zero
							//elements were not written to
		a._avx1 = _mm256_mul_pd(a._avx1,mask);
		__m256d s = _mm256_hadd_pd(a._avx0,a._avx0);
		double new_eri_double =  ((double*)&s)[0] + ((double*)&s)[2];
		s = _mm256_hadd_pd(a._avx1,a._avx1);
		new_eri_double +=  ((double*)&s)[0] + ((double*)&s)[2];
		return new_eri_double;
	}	

template <typename Real>
class VectorSIMD<Real,6>{
public:
		__m256d _avx0;
		__m256d _avx1;
		VectorSIMD(){};
		VectorSIMD(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			_avx1=_mm256_set_pd(0,0,a,a);
		}

		VectorSIMD(Real (&a)[6]){
			_avx0=_mm256_set_pd(a[3],a[2],a[1],a[0]);
			_avx1=_mm256_set_pd(0,0,a[5],a[4]);
		}

		VectorSIMD(Real _0,Real _1,Real _2,Real _3,Real _4,Real _5){
			_avx0=_mm256_set_pd(_3,_2,_1,_0);
			_avx1=_mm256_set_pd( 0,0,_5,_4);
		}

		VectorSIMD(__m256d _0,__m256d _1){
			_avx0= _0;
			_avx1= _1;
		}

		VectorSIMD& operator=(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			_avx1 = _mm256_set_pd(0,0,a,a);
			return *this;
		}

		VectorSIMD& operator += (VectorSIMD<Real,6> a){
			_avx0=  _mm256_add_pd(_avx0,a._avx0);
			_avx1=  _mm256_add_pd(_avx1,a._avx1);
			return *this;
		}

		VectorSIMD& operator -= (VectorSIMD<Real,6> a){
			_avx0=  _mm256_sub_pd(_avx0,a._avx0);
			_avx1=  _mm256_sub_pd(_avx1,a._avx1);
			return *this;
		}

		VectorSIMD operator -() const{
			const static __m256d m1 = _mm256_set_pd(-1.0,-1.0,-1.0,-1.0);
			VectorSIMD<double,6> result;
			result._avx0=  _mm256_mul_pd(this->_avx0,m1);
			result._avx1=  _mm256_mul_pd(this->_avx1,m1);
			return result;
		}

		void load(Real const* a){
			_avx0 =  _mm256_loadu_pd(&a[0]);
			_avx1 =  _mm256_loadu_pd(&a[4]);
		}

		void load_aligned(Real const* a){
			_avx0 =  _mm256_load_pd(&a[0]);
			_avx1 =  _mm256_load_pd(&a[4]);
		}

		void convert(Real *a) const {
			double temp[8];
			_mm256_storeu_pd(&temp[0],_avx0);
			_mm256_storeu_pd(&temp[4],_avx1);
			a[0] = temp[0];
			a[1] = temp[1];
			a[2] = temp[2];
			a[3] = temp[3];
			a[4] = temp[4];
			a[5] = temp[5];
		}

		void convert_aligned(double *a) const {
			double temp[8];
			_mm256_storeu_pd(&temp[0],_avx0);
			_mm256_storeu_pd(&temp[4],_avx1);
			a[0] = temp[0];
			a[1] = temp[1];
			a[2] = temp[2];
			a[3] = temp[3];
			a[4] = temp[4];
			a[5] = temp[5];
		}

	};


	inline VectorSIMD<double,6> operator*(double a, VectorSIMD<double,6> b){
		VectorSIMD<double,6> c;
		__m256d _a =  _mm256_set_pd(a,a,a,a);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		c._avx1=  _mm256_mul_pd(_a, b._avx1);
		return c;
	}


	inline VectorSIMD<double,6> operator*(VectorSIMD<double,6> a ,double b){
		VectorSIMD<double,6> c;
		__m256d _b =  _mm256_set_pd(b,b,b,b);
		c._avx0=  _mm256_mul_pd(a._avx0, _b);
		c._avx1=  _mm256_mul_pd(a._avx1, _b);
		return c;
	}


	inline VectorSIMD<double,6> operator*(int a, VectorSIMD<double,6> b){
		if(a==1){return b;}
		VectorSIMD<double,6> c;
		double q = static_cast<double>(a);
		__m256d _a =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		c._avx1=  _mm256_mul_pd(_a, b._avx1);
		return c;
	}


	inline VectorSIMD<double,6> operator*(VectorSIMD<double,6> a, int b){
		if(b==1){return a;}
		VectorSIMD<double,6> c;
		double q = static_cast<double>(b);
		__m256d _b =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_b, a._avx0);
		c._avx1=  _mm256_mul_pd(_b, a._avx1);
		return c;
	}


	inline VectorSIMD<double,6> operator*(VectorSIMD<double,6> a, VectorSIMD<double,6> b){
		VectorSIMD<double,6> c;
		c._avx0=  _mm256_mul_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_mul_pd(a._avx1, b._avx1);
		return c;
	}


	inline VectorSIMD<double,6> operator+(VectorSIMD<double,6> a, VectorSIMD<double,6> b){
		VectorSIMD<double,6> c;
		c._avx0=  _mm256_add_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_add_pd(a._avx1, b._avx1);
		return c;
	}


	inline VectorSIMD<double,6> operator-(VectorSIMD<double,6> a, VectorSIMD<double,6> b){
		VectorSIMD<double,6> c;
		c._avx0=  _mm256_sub_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_sub_pd(a._avx1, b._avx1);
		return c;
	}


	inline VectorSIMD<double,6> operator/(VectorSIMD<double,6> a, VectorSIMD<double,6> b){
		VectorSIMD<double,6> c;
		c._avx0=  _mm256_div_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_div_pd(a._avx1, b._avx1);
		return c;
	}


	inline VectorSIMD<double,6> exp(VectorSIMD<double,6> a){
		double a_d[6]; a.convert(a_d);
		a_d[0] = std::exp(a_d[0]);
		a_d[1] = std::exp(a_d[1]);
		a_d[2] = std::exp(a_d[2]);
		a_d[3] = std::exp(a_d[3]);
		a_d[4] = std::exp(a_d[4]);
		a_d[5] = std::exp(a_d[5]);
		VectorSIMD<double,6> r(a_d);
		return r;
	}

	inline VectorSIMD<double,6> sqrt(VectorSIMD<double,6> a){
		double a_d[6]; a.convert(a_d);
		a_d[0] = std::sqrt(a_d[0]);
		a_d[1] = std::sqrt(a_d[1]);
		a_d[2] = std::sqrt(a_d[2]);
		a_d[3] = std::sqrt(a_d[3]);
		a_d[4] = std::sqrt(a_d[4]);
		a_d[5] = std::sqrt(a_d[5]);
		VectorSIMD<double,6> r(a_d);
		return r;
	}

	inline std::ostream& operator<<(std::ostream& os, VectorSIMD<double,6> a){ 
		double ad[6];
		a.convert(ad);
		os << "{" << ad[0] << ","  << ad[1] << ","  << ad[2] << ","  << ad[3] << ","  << ad[4] << ","  << ad[5] << "}";
		return os;
		}

	inline double reduce(VectorSIMD<double,6> a){
		__m256d mask = _mm256_set_pd(0,0,1,1); //this ensures that assumed zero
							//elements were not written to
		a._avx1 = _mm256_mul_pd(a._avx1,mask);
		__m256d s = _mm256_hadd_pd(a._avx0,a._avx0);
		double new_eri_double =  ((double*)&s)[0] + ((double*)&s)[2];
		s = _mm256_hadd_pd(a._avx1,a._avx1);
		new_eri_double +=  ((double*)&s)[0] + ((double*)&s)[2];
		return new_eri_double;
	}	

template <typename Real>
class VectorSIMD<Real,7>{
public:
		__m256d _avx0;
		__m256d _avx1;
		VectorSIMD(){};
		VectorSIMD(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			_avx1=_mm256_set_pd(0,a,a,a);
		}

		VectorSIMD(Real (&a)[7]){
			_avx0=_mm256_set_pd(a[3],a[2],a[1],a[0]);
			_avx1=_mm256_set_pd(0,a[6],a[5],a[4]);
		}

		VectorSIMD(Real _0,Real _1,Real _2,Real _3,Real _4,Real _5,Real _6){
			_avx0=_mm256_set_pd(_3,_2,_1,_0);
			_avx1=_mm256_set_pd( 0,_6,_5,_4);
		}

		VectorSIMD(__m256d _0,__m256d _1){
			_avx0= _0;
			_avx1= _1;
		}

		VectorSIMD& operator=(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			_avx1=_mm256_set_pd(a,a,a,a);
			return *this;
		}

		VectorSIMD& operator += (VectorSIMD<Real,7> a){
			_avx0=  _mm256_add_pd(_avx0,a._avx0);
			_avx1=  _mm256_add_pd(_avx1,a._avx1);
			return *this;
		}

		VectorSIMD& operator -= (VectorSIMD<Real,7> a){
			_avx0=  _mm256_sub_pd(_avx0,a._avx0);
			_avx1=  _mm256_sub_pd(_avx1,a._avx1);
			return *this;
		}

		VectorSIMD operator -() const{
			const static __m256d m1 = _mm256_set_pd(-1.0,-1.0,-1.0,-1.0);
			VectorSIMD<double,7> result;
			result._avx0=  _mm256_mul_pd(this->_avx0,m1);
			result._avx1=  _mm256_mul_pd(this->_avx1,m1);
			return result;
		}

		void load(Real const* a){
			_avx0 =  _mm256_loadu_pd(&a[0]);
			_avx1 =  _mm256_loadu_pd(&a[4]);
		}

		void load_aligned(Real const* a){
			_avx0 =  _mm256_load_pd(&a[0]);
			_avx1 =  _mm256_load_pd(&a[4]);
		}

		void convert(Real *a) const {
			_mm256_storeu_pd(&a[0],_avx0);
			_mm256_storeu_pd(&a[4],_avx1);
		}

		void convert_aligned(Real *a) const {
			_mm256_storeu_pd(&a[0],_avx0);
			_mm256_storeu_pd(&a[4],_avx1);
		}

	};


	inline VectorSIMD<double,7> operator*(double a, VectorSIMD<double,7> b){
		VectorSIMD<double,7> c;
		__m256d _a =  _mm256_set_pd(a,a,a,a);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		c._avx1=  _mm256_mul_pd(_a, b._avx1);
		return c;
	}


	inline VectorSIMD<double,7> operator*(VectorSIMD<double,7> a ,double b){
		VectorSIMD<double,7> c;
		__m256d _b =  _mm256_set_pd(b,b,b,b);
		c._avx0=  _mm256_mul_pd(a._avx0, _b);
		c._avx1=  _mm256_mul_pd(a._avx1, _b);
		return c;
	}


	inline VectorSIMD<double,7> operator*(int a, VectorSIMD<double,7> b){
		if(a==1){return b;}
		VectorSIMD<double,7> c;
		double q = static_cast<double>(a);
		__m256d _a =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		c._avx1=  _mm256_mul_pd(_a, b._avx1);
		return c;
	}


	inline VectorSIMD<double,7> operator*(VectorSIMD<double,7> a, int b){
		if(b==1){return a;}
		VectorSIMD<double,7> c;
		double q = static_cast<double>(b);
		__m256d _b =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_b, a._avx0);
		c._avx1=  _mm256_mul_pd(_b, a._avx1);
		return c;
	}


	inline VectorSIMD<double,7> operator*(VectorSIMD<double,7> a, VectorSIMD<double,7> b){
		VectorSIMD<double,7> c;
		c._avx0=  _mm256_mul_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_mul_pd(a._avx1, b._avx1);
		return c;
	}


	inline VectorSIMD<double,7> operator+(VectorSIMD<double,7> a, VectorSIMD<double,7> b){
		VectorSIMD<double,7> c;
		c._avx0=  _mm256_add_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_add_pd(a._avx1, b._avx1);
		return c;
	}


	inline VectorSIMD<double,7> operator-(VectorSIMD<double,7> a, VectorSIMD<double,7> b){
		VectorSIMD<double,7> c;
		c._avx0=  _mm256_sub_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_sub_pd(a._avx1, b._avx1);
		return c;
	}


	inline VectorSIMD<double,7> operator/(VectorSIMD<double,7> a, VectorSIMD<double,7> b){
		VectorSIMD<double,7> c;
		c._avx0=  _mm256_div_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_div_pd(a._avx1, b._avx1);
		return c;
	}


	inline VectorSIMD<double,7> exp(VectorSIMD<double,7> a){
		double a_d[7]; a.convert(a_d);
		a_d[0] = std::exp(a_d[0]);
		a_d[1] = std::exp(a_d[1]);
		a_d[2] = std::exp(a_d[2]);
		a_d[3] = std::exp(a_d[3]);
		a_d[4] = std::exp(a_d[4]);
		a_d[5] = std::exp(a_d[5]);
		a_d[6] = std::exp(a_d[6]);
		VectorSIMD<double,7> r(a_d);
		return r;
	}

	inline VectorSIMD<double,7> sqrt(VectorSIMD<double,7> a){
		double a_d[7]; a.convert(a_d);
		a_d[0] = std::sqrt(a_d[0]);
		a_d[1] = std::sqrt(a_d[1]);
		a_d[2] = std::sqrt(a_d[2]);
		a_d[3] = std::sqrt(a_d[3]);
		a_d[4] = std::sqrt(a_d[4]);
		a_d[5] = std::sqrt(a_d[5]);
		a_d[6] = std::sqrt(a_d[6]);
		VectorSIMD<double,7> r(a_d);
		return r;
	}

	inline std::ostream& operator<<(std::ostream& os, VectorSIMD<double,7> a){ 
		double ad[7];
		a.convert(ad);
		os << "{" << ad[0] << ","  << ad[1] << ","  << ad[2] << ","  << ad[3] << ","  << ad[4] << ","  << ad[5] << ","  << ad[6] << "}";
		return os;
		}

	inline double reduce(VectorSIMD<double,7> a){
		__m256d mask = _mm256_set_pd(0,1,1,1); //this ensures that assumed zero
							//elements were not written to
		a._avx1 = _mm256_mul_pd(a._avx1,mask);
		__m256d s = _mm256_hadd_pd(a._avx0,a._avx0);
		double new_eri_double =  ((double*)&s)[0] + ((double*)&s)[2];
		s = _mm256_hadd_pd(a._avx1,a._avx1);
		new_eri_double +=  ((double*)&s)[0] + ((double*)&s)[2];
		return new_eri_double;
	}	


template <typename Real>
class VectorSIMD<Real,8>{
public:
		__m256d _avx0;
		__m256d _avx1;
		VectorSIMD(){};
		VectorSIMD(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			_avx1=_mm256_set_pd(a,a,a,a);
		}

		VectorSIMD(Real (&a)[8]){
			_avx0=_mm256_set_pd(a[3],a[2],a[1],a[0]);
			_avx1=_mm256_set_pd(a[7],a[6],a[5],a[4]);
		}

		VectorSIMD(Real _0,Real _1,Real _2,Real _3,Real _4,Real _5,Real _6,Real _7){
			_avx0=_mm256_set_pd(_3,_2,_1,_0);
			_avx1=_mm256_set_pd(_7,_6,_5,_4);
		}

		VectorSIMD(__m256d _0,__m256d _1){
			_avx0= _0;
			_avx1= _1;
		}

		VectorSIMD& operator=(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			_avx1=_mm256_set_pd(a,a,a,a);
			return *this;
		}

		VectorSIMD& operator += (VectorSIMD<Real,8> a){
			_avx0=  _mm256_add_pd(_avx0,a._avx0);
			_avx1=  _mm256_add_pd(_avx1,a._avx1);
			return *this;
		}

		VectorSIMD& operator -= (VectorSIMD<Real,8> a){
			_avx0=  _mm256_sub_pd(_avx0,a._avx0);
			_avx1=  _mm256_sub_pd(_avx1,a._avx1);
			return *this;
		}

		VectorSIMD operator -() const{
			const static __m256d m1 = _mm256_set_pd(-1.0,-1.0,-1.0,-1.0);
			VectorSIMD<double,8> result;
			result._avx0=  _mm256_mul_pd(this->_avx0,m1);
			result._avx1=  _mm256_mul_pd(this->_avx1,m1);
			return result;
		}

		void load(Real const* a){
			_avx0 =  _mm256_loadu_pd(&a[0]);
			_avx1 =  _mm256_loadu_pd(&a[4]);
		}

		void load_aligned(Real const* a){
			_avx0 =  _mm256_load_pd(&a[0]);
			_avx1 =  _mm256_load_pd(&a[4]);
		}

		void convert(Real *a) const {
			_mm256_storeu_pd(&a[0],_avx0);
			_mm256_storeu_pd(&a[4],_avx1);
		}

		void convert_aligned(Real *a) const {
			_mm256_storeu_pd(&a[0],_avx0);
			_mm256_storeu_pd(&a[4],_avx1);
		}

	};


	inline VectorSIMD<double,8> operator*(double a, VectorSIMD<double,8> b){
		VectorSIMD<double,8> c;
		__m256d _a =  _mm256_set_pd(a,a,a,a);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		c._avx1=  _mm256_mul_pd(_a, b._avx1);
		return c;
	}


	inline VectorSIMD<double,8> operator*(VectorSIMD<double,8> a ,double b){
		VectorSIMD<double,8> c;
		__m256d _b =  _mm256_set_pd(b,b,b,b);
		c._avx0=  _mm256_mul_pd(a._avx0, _b);
		c._avx1=  _mm256_mul_pd(a._avx1, _b);
		return c;
	}


	inline VectorSIMD<double,8> operator*(int a, VectorSIMD<double,8> b){
		if(a==1){return b;}
		VectorSIMD<double,8> c;
		double q = static_cast<double>(a);
		__m256d _a =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		c._avx1=  _mm256_mul_pd(_a, b._avx1);
		return c;
	}


	inline VectorSIMD<double,8> operator*(VectorSIMD<double,8> a, int b){
		if(b==1){return a;}
		VectorSIMD<double,8> c;
		double q = static_cast<double>(b);
		__m256d _b =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_b, a._avx0);
		c._avx1=  _mm256_mul_pd(_b, a._avx1);
		return c;
	}


	inline VectorSIMD<double,8> operator*(VectorSIMD<double,8> a, VectorSIMD<double,8> b){
		VectorSIMD<double,8> c;
		c._avx0=  _mm256_mul_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_mul_pd(a._avx1, b._avx1);
		return c;
	}


	inline VectorSIMD<double,8> operator+(VectorSIMD<double,8> a, VectorSIMD<double,8> b){
		VectorSIMD<double,8> c;
		c._avx0=  _mm256_add_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_add_pd(a._avx1, b._avx1);
		return c;
	}


	inline VectorSIMD<double,8> operator-(VectorSIMD<double,8> a, VectorSIMD<double,8> b){
		VectorSIMD<double,8> c;
		c._avx0=  _mm256_sub_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_sub_pd(a._avx1, b._avx1);
		return c;
	}


	inline VectorSIMD<double,8> operator/(VectorSIMD<double,8> a, VectorSIMD<double,8> b){
		VectorSIMD<double,8> c;
		c._avx0=  _mm256_div_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_div_pd(a._avx1, b._avx1);
		return c;
	}


	inline VectorSIMD<double,8> exp(VectorSIMD<double,8> a){
		double a_d[8]; a.convert(a_d);
		a_d[0] = std::exp(a_d[0]);
		a_d[1] = std::exp(a_d[1]);
		a_d[2] = std::exp(a_d[2]);
		a_d[3] = std::exp(a_d[3]);
		a_d[4] = std::exp(a_d[4]);
		a_d[5] = std::exp(a_d[5]);
		a_d[6] = std::exp(a_d[6]);
		a_d[7] = std::exp(a_d[7]);
		VectorSIMD<double,8> r(a_d);
		return r;
	}

	inline VectorSIMD<double,8> sqrt(VectorSIMD<double,8> a){
		double a_d[8]; a.convert(a_d);
		a_d[0] = std::sqrt(a_d[0]);
		a_d[1] = std::sqrt(a_d[1]);
		a_d[2] = std::sqrt(a_d[2]);
		a_d[3] = std::sqrt(a_d[3]);
		a_d[4] = std::sqrt(a_d[4]);
		a_d[5] = std::sqrt(a_d[5]);
		a_d[6] = std::sqrt(a_d[6]);
		a_d[7] = std::sqrt(a_d[7]);
		VectorSIMD<double,8> r(a_d);
		return r;
	}

	inline std::ostream& operator<<(std::ostream& os, VectorSIMD<double,8> a){ 
		double ad[8];
		a.convert(ad);
		os << "{" << ad[0] << ","  << ad[1] << ","  << ad[2] << ","  << ad[3] << ","  << ad[4] << ","  << ad[5] << ","  << ad[6] << ","  << ad[7] << "}";
		return os;
		}

	inline double reduce(VectorSIMD<double,8> a){
		__m256d s = _mm256_hadd_pd(a._avx0,a._avx0);
		double new_eri_double =  ((double*)&s)[0] + ((double*)&s)[2];
		s = _mm256_hadd_pd(a._avx1,a._avx1);
		new_eri_double +=  ((double*)&s)[0] + ((double*)&s)[2];
		return new_eri_double;
	}	

template <typename Real>
class VectorSIMD<Real,9>{
public:
		__m256d _avx0;
		__m256d _avx1;
		__m256d _avx2;
		VectorSIMD(){};
		VectorSIMD(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			_avx1=_mm256_set_pd(a,a,a,a);
			_avx2=_mm256_set_pd(0,0,0,a);
		}

		VectorSIMD(Real (&a)[9]){
			_avx0=_mm256_set_pd(a[3],a[2],a[1],a[0]);
			_avx1=_mm256_set_pd(a[7],a[6],a[5],a[4]);
			_avx2=_mm256_set_pd(0,0,0,a[8]);
		}

		VectorSIMD(Real _0,Real _1,Real _2,Real _3,Real _4,Real _5,Real _6,Real _7,Real _8){
			_avx0=_mm256_set_pd(_3,_2,_1,_0);
			_avx1=_mm256_set_pd(_7,_6,_5,_4);
			_avx2=_mm256_set_pd( 0,0,0,_8);
		}

		VectorSIMD(__m256d _0,__m256d _1,__m256d _2){
			_avx0= _0;
			_avx1= _1;
			_avx2= _2;
		}

		VectorSIMD& operator=(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			_avx1=_mm256_set_pd(a,a,a,a);
			_avx2 = _mm256_set_pd(0,0,0,a);
			return *this;
		}

		VectorSIMD& operator += (VectorSIMD<Real,9> a){
			_avx0=  _mm256_add_pd(_avx0,a._avx0);
			_avx1=  _mm256_add_pd(_avx1,a._avx1);
			_avx2=  _mm256_add_pd(_avx2,a._avx2);
			return *this;
		}

		VectorSIMD& operator -= (VectorSIMD<Real,9> a){
			_avx0=  _mm256_sub_pd(_avx0,a._avx0);
			_avx1=  _mm256_sub_pd(_avx1,a._avx1);
			_avx2=  _mm256_sub_pd(_avx2,a._avx2);
			return *this;
		}

		VectorSIMD operator -() const{
			const static __m256d m1 = _mm256_set_pd(-1.0,-1.0,-1.0,-1.0);
			VectorSIMD<double,9> result;
			result._avx0=  _mm256_mul_pd(this->_avx0,m1);
			result._avx1=  _mm256_mul_pd(this->_avx1,m1);
			result._avx2=  _mm256_mul_pd(this->_avx2,m1);
			return result;
		}

		void load(Real const* a){
			_avx0 =  _mm256_loadu_pd(&a[0]);
			_avx1 =  _mm256_loadu_pd(&a[4]);
			_avx2 =  _mm256_loadu_pd(&a[8]);
		}

		void load_aligned(Real const* a){
			_avx0 =  _mm256_load_pd(&a[0]);
			_avx1 =  _mm256_load_pd(&a[4]);
			_avx2 =  _mm256_load_pd(&a[8]);
		}

		void convert(Real *a) const {
			double temp[12];
			_mm256_storeu_pd(&temp[0],_avx0);
			_mm256_storeu_pd(&temp[4],_avx1);
			_mm256_storeu_pd(&temp[8],_avx2);
			a[0] = temp[0];
			a[1] = temp[1];
			a[2] = temp[2];
			a[3] = temp[3];
			a[4] = temp[4];
			a[5] = temp[5];
			a[6] = temp[6];
			a[7] = temp[7];
			a[8] = temp[8];
		}

		void convert_aligned(double *a) const {
			double temp[12];
			_mm256_storeu_pd(&temp[0],_avx0);
			_mm256_storeu_pd(&temp[4],_avx1);
			_mm256_storeu_pd(&temp[8],_avx2);
			a[0] = temp[0];
			a[1] = temp[1];
			a[2] = temp[2];
			a[3] = temp[3];
			a[4] = temp[4];
			a[5] = temp[5];
			a[6] = temp[6];
			a[7] = temp[7];
			a[8] = temp[8];
		}

	};


	inline VectorSIMD<double,9> operator*(double a, VectorSIMD<double,9> b){
		VectorSIMD<double,9> c;
		__m256d _a =  _mm256_set_pd(a,a,a,a);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		c._avx1=  _mm256_mul_pd(_a, b._avx1);
		c._avx2=  _mm256_mul_pd(_a, b._avx2);
		return c;
	}


	inline VectorSIMD<double,9> operator*(VectorSIMD<double,9> a ,double b){
		VectorSIMD<double,9> c;
		__m256d _b =  _mm256_set_pd(b,b,b,b);
		c._avx0=  _mm256_mul_pd(a._avx0, _b);
		c._avx1=  _mm256_mul_pd(a._avx1, _b);
		c._avx2=  _mm256_mul_pd(a._avx2, _b);
		return c;
	}


	inline VectorSIMD<double,9> operator*(int a, VectorSIMD<double,9> b){
		if(a==1){return b;}
		VectorSIMD<double,9> c;
		double q = static_cast<double>(a);
		__m256d _a =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		c._avx1=  _mm256_mul_pd(_a, b._avx1);
		c._avx2=  _mm256_mul_pd(_a, b._avx2);
		return c;
	}


	inline VectorSIMD<double,9> operator*(VectorSIMD<double,9> a, int b){
		if(b==1){return a;}
		VectorSIMD<double,9> c;
		double q = static_cast<double>(b);
		__m256d _b =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_b, a._avx0);
		c._avx1=  _mm256_mul_pd(_b, a._avx1);
		c._avx2=  _mm256_mul_pd(_b, a._avx2);
		return c;
	}


	inline VectorSIMD<double,9> operator*(VectorSIMD<double,9> a, VectorSIMD<double,9> b){
		VectorSIMD<double,9> c;
		c._avx0=  _mm256_mul_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_mul_pd(a._avx1, b._avx1);
		c._avx2=  _mm256_mul_pd(a._avx2, b._avx2);
		return c;
	}


	inline VectorSIMD<double,9> operator+(VectorSIMD<double,9> a, VectorSIMD<double,9> b){
		VectorSIMD<double,9> c;
		c._avx0=  _mm256_add_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_add_pd(a._avx1, b._avx1);
		c._avx2=  _mm256_add_pd(a._avx2, b._avx2);
		return c;
	}


	inline VectorSIMD<double,9> operator-(VectorSIMD<double,9> a, VectorSIMD<double,9> b){
		VectorSIMD<double,9> c;
		c._avx0=  _mm256_sub_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_sub_pd(a._avx1, b._avx1);
		c._avx2=  _mm256_sub_pd(a._avx2, b._avx2);
		return c;
	}


	inline VectorSIMD<double,9> operator/(VectorSIMD<double,9> a, VectorSIMD<double,9> b){
		VectorSIMD<double,9> c;
		c._avx0=  _mm256_div_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_div_pd(a._avx1, b._avx1);
		c._avx2=  _mm256_div_pd(a._avx2, b._avx2);
		return c;
	}


	inline VectorSIMD<double,9> exp(VectorSIMD<double,9> a){
		double a_d[9]; a.convert(a_d);
		a_d[0] = std::exp(a_d[0]);
		a_d[1] = std::exp(a_d[1]);
		a_d[2] = std::exp(a_d[2]);
		a_d[3] = std::exp(a_d[3]);
		a_d[4] = std::exp(a_d[4]);
		a_d[5] = std::exp(a_d[5]);
		a_d[6] = std::exp(a_d[6]);
		a_d[7] = std::exp(a_d[7]);
		a_d[8] = std::exp(a_d[8]);
		VectorSIMD<double,9> r(a_d);
		return r;
	}

	inline VectorSIMD<double,9> sqrt(VectorSIMD<double,9> a){
		double a_d[9]; a.convert(a_d);
		a_d[0] = std::sqrt(a_d[0]);
		a_d[1] = std::sqrt(a_d[1]);
		a_d[2] = std::sqrt(a_d[2]);
		a_d[3] = std::sqrt(a_d[3]);
		a_d[4] = std::sqrt(a_d[4]);
		a_d[5] = std::sqrt(a_d[5]);
		a_d[6] = std::sqrt(a_d[6]);
		a_d[7] = std::sqrt(a_d[7]);
		a_d[8] = std::sqrt(a_d[8]);
		VectorSIMD<double,9> r(a_d);
		return r;
	}

	inline std::ostream& operator<<(std::ostream& os, VectorSIMD<double,9> a){ 
		double ad[9];
		a.convert(ad);
		os << "{" << ad[0] << ","  << ad[1] << ","  << ad[2] << ","  << ad[3] << ","  << ad[4] << ","  << ad[5] << ","  << ad[6] << ","  << ad[7] << ","  << ad[8] << "}";
		return os;
		}


	inline double reduce(VectorSIMD<double,9> a){
		__m256d mask = _mm256_set_pd(0,0,0,1); //this ensures that assumed zero
							//elements were not written to
		a._avx2 = _mm256_mul_pd(a._avx2,mask);
		__m256d s = _mm256_hadd_pd(a._avx0,a._avx0);
		double new_eri_double =  ((double*)&s)[0] + ((double*)&s)[2];
		s = _mm256_hadd_pd(a._avx1,a._avx1);
		new_eri_double +=  ((double*)&s)[0] + ((double*)&s)[2];
		s = _mm256_hadd_pd(a._avx2,a._avx2);
		new_eri_double +=  ((double*)&s)[0] + ((double*)&s)[2];
		return new_eri_double;
	}	



template <typename Real>
class VectorSIMD<Real,10>{
public:
		__m256d _avx0;
		__m256d _avx1;
		__m256d _avx2;
		VectorSIMD(){};
		VectorSIMD(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			_avx1=_mm256_set_pd(a,a,a,a);
			_avx2=_mm256_set_pd(0,0,a,a);
		}

		VectorSIMD(Real (&a)[10]){
			_avx0=_mm256_set_pd(a[3],a[2],a[1],a[0]);
			_avx1=_mm256_set_pd(a[7],a[6],a[5],a[4]);
			_avx2=_mm256_set_pd(0,0,a[9],a[8]);
		}

		VectorSIMD(Real _0,Real _1,Real _2,Real _3,Real _4,Real _5,Real _6,Real _7,Real _8,Real _9){
			_avx0=_mm256_set_pd(_3,_2,_1,_0);
			_avx1=_mm256_set_pd(_7,_6,_5,_4);
			_avx2=_mm256_set_pd( 0,0,_9,_8);
		}

		VectorSIMD(__m256d _0,__m256d _1,__m256d _2){
			_avx0= _0;
			_avx1= _1;
			_avx2= _2;
		}

		VectorSIMD& operator=(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			_avx1=_mm256_set_pd(a,a,a,a);
			_avx2 = _mm256_set_pd(0,0,a,a);
			return *this;
		}

		VectorSIMD& operator += (VectorSIMD<Real,10> a){
			_avx0=  _mm256_add_pd(_avx0,a._avx0);
			_avx1=  _mm256_add_pd(_avx1,a._avx1);
			_avx2=  _mm256_add_pd(_avx2,a._avx2);
			return *this;
		}

		VectorSIMD& operator -= (VectorSIMD<Real,10> a){
			_avx0=  _mm256_sub_pd(_avx0,a._avx0);
			_avx1=  _mm256_sub_pd(_avx1,a._avx1);
			_avx2=  _mm256_sub_pd(_avx2,a._avx2);
			return *this;
		}

		VectorSIMD operator -() const{
			const static __m256d m1 = _mm256_set_pd(-1.0,-1.0,-1.0,-1.0);
			VectorSIMD<double,10> result;
			result._avx0=  _mm256_mul_pd(this->_avx0,m1);
			result._avx1=  _mm256_mul_pd(this->_avx1,m1);
			result._avx2=  _mm256_mul_pd(this->_avx2,m1);
			return result;
		}

		void load(Real const* a){
			_avx0 =  _mm256_loadu_pd(&a[0]);
			_avx1 =  _mm256_loadu_pd(&a[4]);
			_avx2 =  _mm256_loadu_pd(&a[8]);
		}

		void load_aligned(Real const* a){
			_avx0 =  _mm256_load_pd(&a[0]);
			_avx1 =  _mm256_load_pd(&a[4]);
			_avx2 =  _mm256_load_pd(&a[8]);
		}

		void convert(Real *a) const {
			double temp[12];
			_mm256_storeu_pd(&temp[0],_avx0);
			_mm256_storeu_pd(&temp[4],_avx1);
			_mm256_storeu_pd(&temp[8],_avx2);
			a[0] = temp[0];
			a[1] = temp[1];
			a[2] = temp[2];
			a[3] = temp[3];
			a[4] = temp[4];
			a[5] = temp[5];
			a[6] = temp[6];
			a[7] = temp[7];
			a[8] = temp[8];
			a[9] = temp[9];
		}

		void convert_aligned(double *a) const {
			double temp[12];
			_mm256_storeu_pd(&temp[0],_avx0);
			_mm256_storeu_pd(&temp[4],_avx1);
			_mm256_storeu_pd(&temp[8],_avx2);
			a[0] = temp[0];
			a[1] = temp[1];
			a[2] = temp[2];
			a[3] = temp[3];
			a[4] = temp[4];
			a[5] = temp[5];
			a[6] = temp[6];
			a[7] = temp[7];
			a[8] = temp[8];
			a[9] = temp[9];
		}

	};


	inline VectorSIMD<double,10> operator*(double a, VectorSIMD<double,10> b){
		VectorSIMD<double,10> c;
		__m256d _a =  _mm256_set_pd(a,a,a,a);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		c._avx1=  _mm256_mul_pd(_a, b._avx1);
		c._avx2=  _mm256_mul_pd(_a, b._avx2);
		return c;
	}


	inline VectorSIMD<double,10> operator*(VectorSIMD<double,10> a ,double b){
		VectorSIMD<double,10> c;
		__m256d _b =  _mm256_set_pd(b,b,b,b);
		c._avx0=  _mm256_mul_pd(a._avx0, _b);
		c._avx1=  _mm256_mul_pd(a._avx1, _b);
		c._avx2=  _mm256_mul_pd(a._avx2, _b);
		return c;
	}


	inline VectorSIMD<double,10> operator*(int a, VectorSIMD<double,10> b){
		if(a==1){return b;}
		VectorSIMD<double,10> c;
		double q = static_cast<double>(a);
		__m256d _a =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		c._avx1=  _mm256_mul_pd(_a, b._avx1);
		c._avx2=  _mm256_mul_pd(_a, b._avx2);
		return c;
	}


	inline VectorSIMD<double,10> operator*(VectorSIMD<double,10> a, int b){
		if(b==1){return a;}
		VectorSIMD<double,10> c;
		double q = static_cast<double>(b);
		__m256d _b =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_b, a._avx0);
		c._avx1=  _mm256_mul_pd(_b, a._avx1);
		c._avx2=  _mm256_mul_pd(_b, a._avx2);
		return c;
	}


	inline VectorSIMD<double,10> operator*(VectorSIMD<double,10> a, VectorSIMD<double,10> b){
		VectorSIMD<double,10> c;
		c._avx0=  _mm256_mul_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_mul_pd(a._avx1, b._avx1);
		c._avx2=  _mm256_mul_pd(a._avx2, b._avx2);
		return c;
	}


	inline VectorSIMD<double,10> operator+(VectorSIMD<double,10> a, VectorSIMD<double,10> b){
		VectorSIMD<double,10> c;
		c._avx0=  _mm256_add_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_add_pd(a._avx1, b._avx1);
		c._avx2=  _mm256_add_pd(a._avx2, b._avx2);
		return c;
	}


	inline VectorSIMD<double,10> operator-(VectorSIMD<double,10> a, VectorSIMD<double,10> b){
		VectorSIMD<double,10> c;
		c._avx0=  _mm256_sub_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_sub_pd(a._avx1, b._avx1);
		c._avx2=  _mm256_sub_pd(a._avx2, b._avx2);
		return c;
	}


	inline VectorSIMD<double,10> operator/(VectorSIMD<double,10> a, VectorSIMD<double,10> b){
		VectorSIMD<double,10> c;
		c._avx0=  _mm256_div_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_div_pd(a._avx1, b._avx1);
		c._avx2=  _mm256_div_pd(a._avx2, b._avx2);
		return c;
	}


	inline VectorSIMD<double,10> exp(VectorSIMD<double,10> a){
		double a_d[10]; a.convert(a_d);
		a_d[0] = std::exp(a_d[0]);
		a_d[1] = std::exp(a_d[1]);
		a_d[2] = std::exp(a_d[2]);
		a_d[3] = std::exp(a_d[3]);
		a_d[4] = std::exp(a_d[4]);
		a_d[5] = std::exp(a_d[5]);
		a_d[6] = std::exp(a_d[6]);
		a_d[7] = std::exp(a_d[7]);
		a_d[8] = std::exp(a_d[8]);
		a_d[9] = std::exp(a_d[9]);
		VectorSIMD<double,10> r(a_d);
		return r;
	}

	inline VectorSIMD<double,10> sqrt(VectorSIMD<double,10> a){
		double a_d[10]; a.convert(a_d);
		a_d[0] = std::sqrt(a_d[0]);
		a_d[1] = std::sqrt(a_d[1]);
		a_d[2] = std::sqrt(a_d[2]);
		a_d[3] = std::sqrt(a_d[3]);
		a_d[4] = std::sqrt(a_d[4]);
		a_d[5] = std::sqrt(a_d[5]);
		a_d[6] = std::sqrt(a_d[6]);
		a_d[7] = std::sqrt(a_d[7]);
		a_d[8] = std::sqrt(a_d[8]);
		a_d[9] = std::sqrt(a_d[9]);
		VectorSIMD<double,10> r(a_d);
		return r;
	}

	inline std::ostream& operator<<(std::ostream& os, VectorSIMD<double,10> a){ 
		double ad[10];
		a.convert(ad);
		os << "{" << ad[0] << ","  << ad[1] << ","  << ad[2] << ","  << ad[3] << ","  << ad[4] << ","  << ad[5] << ","  << ad[6] << ","  << ad[7] << ","  << ad[8] << ","  << ad[9] << "}";
		return os;
		}

	inline double reduce(VectorSIMD<double,10> a){
		__m256d mask = _mm256_set_pd(0,0,1,1); //this ensures that assumed zero
							//elements were not written to
		a._avx2 = _mm256_mul_pd(a._avx2,mask);
		__m256d s = _mm256_hadd_pd(a._avx0,a._avx0);
		double new_eri_double =  ((double*)&s)[0] + ((double*)&s)[2];
		s = _mm256_hadd_pd(a._avx1,a._avx1);
		new_eri_double +=  ((double*)&s)[0] + ((double*)&s)[2];
		s = _mm256_hadd_pd(a._avx2,a._avx2);
		new_eri_double +=  ((double*)&s)[0] + ((double*)&s)[2];
		return new_eri_double;
	}	

template <typename Real>
class VectorSIMD<Real,11>{
public:
		__m256d _avx0;
		__m256d _avx1;
		__m256d _avx2;
		VectorSIMD(){};
		VectorSIMD(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			_avx1=_mm256_set_pd(a,a,a,a);
			_avx2=_mm256_set_pd(0,a,a,a);
		}

		VectorSIMD(Real (&a)[11]){
			_avx0=_mm256_set_pd(a[3],a[2],a[1],a[0]);
			_avx1=_mm256_set_pd(a[7],a[6],a[5],a[4]);
			_avx2=_mm256_set_pd(0,a[10],a[9],a[8]);
		}

		VectorSIMD(Real _0,Real _1,Real _2,Real _3,Real _4,Real _5,Real _6,Real _7,Real _8,Real _9,Real _10){
			_avx0=_mm256_set_pd(_3,_2,_1,_0);
			_avx1=_mm256_set_pd(_7,_6,_5,_4);
			_avx2=_mm256_set_pd( 0,_10,_9,_8);
		}

		VectorSIMD(__m256d _0,__m256d _1,__m256d _2){
			_avx0= _0;
			_avx1= _1;
			_avx2= _2;
		}

		VectorSIMD& operator=(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			_avx1=_mm256_set_pd(a,a,a,a);
			_avx2=_mm256_set_pd(a,a,a,a);
			return *this;
		}

		VectorSIMD& operator += (VectorSIMD<Real,11> a){
			_avx0=  _mm256_add_pd(_avx0,a._avx0);
			_avx1=  _mm256_add_pd(_avx1,a._avx1);
			_avx2=  _mm256_add_pd(_avx2,a._avx2);
			return *this;
		}

		VectorSIMD& operator -= (VectorSIMD<Real,11> a){
			_avx0=  _mm256_sub_pd(_avx0,a._avx0);
			_avx1=  _mm256_sub_pd(_avx1,a._avx1);
			_avx2=  _mm256_sub_pd(_avx2,a._avx2);
			return *this;
		}

		VectorSIMD operator -() const{
			const static __m256d m1 = _mm256_set_pd(-1.0,-1.0,-1.0,-1.0);
			VectorSIMD<double,11> result;
			result._avx0=  _mm256_mul_pd(this->_avx0,m1);
			result._avx1=  _mm256_mul_pd(this->_avx1,m1);
			result._avx2=  _mm256_mul_pd(this->_avx2,m1);
			return result;
		}

		void load(Real const* a){
			_avx0 =  _mm256_loadu_pd(&a[0]);
			_avx1 =  _mm256_loadu_pd(&a[4]);
			_avx2 =  _mm256_loadu_pd(&a[8]);
		}

		void load_aligned(Real const* a){
			_avx0 =  _mm256_load_pd(&a[0]);
			_avx1 =  _mm256_load_pd(&a[4]);
			_avx2 =  _mm256_load_pd(&a[8]);
		}

		void convert(Real *a) const {
			_mm256_storeu_pd(&a[0],_avx0);
			_mm256_storeu_pd(&a[4],_avx1);
			_mm256_storeu_pd(&a[8],_avx2);
		}

		void convert_aligned(Real *a) const {
			_mm256_storeu_pd(&a[0],_avx0);
			_mm256_storeu_pd(&a[4],_avx1);
			_mm256_storeu_pd(&a[8],_avx2);
		}

	};


	inline VectorSIMD<double,11> operator*(double a, VectorSIMD<double,11> b){
		VectorSIMD<double,11> c;
		__m256d _a =  _mm256_set_pd(a,a,a,a);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		c._avx1=  _mm256_mul_pd(_a, b._avx1);
		c._avx2=  _mm256_mul_pd(_a, b._avx2);
		return c;
	}


	inline VectorSIMD<double,11> operator*(VectorSIMD<double,11> a ,double b){
		VectorSIMD<double,11> c;
		__m256d _b =  _mm256_set_pd(b,b,b,b);
		c._avx0=  _mm256_mul_pd(a._avx0, _b);
		c._avx1=  _mm256_mul_pd(a._avx1, _b);
		c._avx2=  _mm256_mul_pd(a._avx2, _b);
		return c;
	}


	inline VectorSIMD<double,11> operator*(int a, VectorSIMD<double,11> b){
		if(a==1){return b;}
		VectorSIMD<double,11> c;
		double q = static_cast<double>(a);
		__m256d _a =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		c._avx1=  _mm256_mul_pd(_a, b._avx1);
		c._avx2=  _mm256_mul_pd(_a, b._avx2);
		return c;
	}


	inline VectorSIMD<double,11> operator*(VectorSIMD<double,11> a, int b){
		if(b==1){return a;}
		VectorSIMD<double,11> c;
		double q = static_cast<double>(b);
		__m256d _b =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_b, a._avx0);
		c._avx1=  _mm256_mul_pd(_b, a._avx1);
		c._avx2=  _mm256_mul_pd(_b, a._avx2);
		return c;
	}


	inline VectorSIMD<double,11> operator*(VectorSIMD<double,11> a, VectorSIMD<double,11> b){
		VectorSIMD<double,11> c;
		c._avx0=  _mm256_mul_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_mul_pd(a._avx1, b._avx1);
		c._avx2=  _mm256_mul_pd(a._avx2, b._avx2);
		return c;
	}


	inline VectorSIMD<double,11> operator+(VectorSIMD<double,11> a, VectorSIMD<double,11> b){
		VectorSIMD<double,11> c;
		c._avx0=  _mm256_add_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_add_pd(a._avx1, b._avx1);
		c._avx2=  _mm256_add_pd(a._avx2, b._avx2);
		return c;
	}


	inline VectorSIMD<double,11> operator-(VectorSIMD<double,11> a, VectorSIMD<double,11> b){
		VectorSIMD<double,11> c;
		c._avx0=  _mm256_sub_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_sub_pd(a._avx1, b._avx1);
		c._avx2=  _mm256_sub_pd(a._avx2, b._avx2);
		return c;
	}


	inline VectorSIMD<double,11> operator/(VectorSIMD<double,11> a, VectorSIMD<double,11> b){
		VectorSIMD<double,11> c;
		c._avx0=  _mm256_div_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_div_pd(a._avx1, b._avx1);
		c._avx2=  _mm256_div_pd(a._avx2, b._avx2);
		return c;
	}


	inline VectorSIMD<double,11> exp(VectorSIMD<double,11> a){
		double a_d[11]; a.convert(a_d);
		a_d[0] = std::exp(a_d[0]);
		a_d[1] = std::exp(a_d[1]);
		a_d[2] = std::exp(a_d[2]);
		a_d[3] = std::exp(a_d[3]);
		a_d[4] = std::exp(a_d[4]);
		a_d[5] = std::exp(a_d[5]);
		a_d[6] = std::exp(a_d[6]);
		a_d[7] = std::exp(a_d[7]);
		a_d[8] = std::exp(a_d[8]);
		a_d[9] = std::exp(a_d[9]);
		a_d[10] = std::exp(a_d[10]);
		VectorSIMD<double,11> r(a_d);
		return r;
	}

	inline VectorSIMD<double,11> sqrt(VectorSIMD<double,11> a){
		double a_d[11]; a.convert(a_d);
		a_d[0] = std::sqrt(a_d[0]);
		a_d[1] = std::sqrt(a_d[1]);
		a_d[2] = std::sqrt(a_d[2]);
		a_d[3] = std::sqrt(a_d[3]);
		a_d[4] = std::sqrt(a_d[4]);
		a_d[5] = std::sqrt(a_d[5]);
		a_d[6] = std::sqrt(a_d[6]);
		a_d[7] = std::sqrt(a_d[7]);
		a_d[8] = std::sqrt(a_d[8]);
		a_d[9] = std::sqrt(a_d[9]);
		a_d[10] = std::sqrt(a_d[10]);
		VectorSIMD<double,11> r(a_d);
		return r;
	}

	inline std::ostream& operator<<(std::ostream& os, VectorSIMD<double,11> a){ 
		double ad[11];
		a.convert(ad);
		os << "{" << ad[0] << ","  << ad[1] << ","  << ad[2] << ","  << ad[3] << ","  << ad[4] << ","  << ad[5] << ","  << ad[6] << ","  << ad[7] << ","  << ad[8] << ","  << ad[9] << ","  << ad[10] << "}";
		return os;
		}

	inline double reduce(VectorSIMD<double,11> a){
		__m256d mask = _mm256_set_pd(0,1,1,1); //this ensures that assumed zero
							//elements were not written to
		a._avx2 = _mm256_mul_pd(a._avx2,mask);
		__m256d s = _mm256_hadd_pd(a._avx0,a._avx0);
		double new_eri_double =  ((double*)&s)[0] + ((double*)&s)[2];
		s = _mm256_hadd_pd(a._avx1,a._avx1);
		new_eri_double +=  ((double*)&s)[0] + ((double*)&s)[2];
		s = _mm256_hadd_pd(a._avx2,a._avx2);
		new_eri_double +=  ((double*)&s)[0] + ((double*)&s)[2];
		return new_eri_double;
	}	


template <typename Real>
class VectorSIMD<Real,12>{
public:
		__m256d _avx0;
		__m256d _avx1;
		__m256d _avx2;
		VectorSIMD(){};
		VectorSIMD(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			_avx1=_mm256_set_pd(a,a,a,a);
			_avx2=_mm256_set_pd(a,a,a,a);
		}

		VectorSIMD(Real (&a)[12]){
			_avx0=_mm256_set_pd(a[3],a[2],a[1],a[0]);
			_avx1=_mm256_set_pd(a[7],a[6],a[5],a[4]);
			_avx2=_mm256_set_pd(a[11],a[10],a[9],a[8]);
		}

		VectorSIMD(Real _0,Real _1,Real _2,Real _3,Real _4,Real _5,Real _6,Real _7,Real _8,Real _9,Real _10,Real _11){
			_avx0=_mm256_set_pd(_3,_2,_1,_0);
			_avx1=_mm256_set_pd(_7,_6,_5,_4);
			_avx2=_mm256_set_pd(_11,_10,_9,_8);
		}

		VectorSIMD(__m256d _0,__m256d _1,__m256d _2){
			_avx0= _0;
			_avx1= _1;
			_avx2= _2;
		}

		VectorSIMD& operator=(Real a){
			_avx0=_mm256_set_pd(a,a,a,a);
			_avx1=_mm256_set_pd(a,a,a,a);
			_avx2=_mm256_set_pd(a,a,a,a);
			return *this;
		}

		VectorSIMD& operator += (VectorSIMD<Real,12> a){
			_avx0=  _mm256_add_pd(_avx0,a._avx0);
			_avx1=  _mm256_add_pd(_avx1,a._avx1);
			_avx2=  _mm256_add_pd(_avx2,a._avx2);
			return *this;
		}

		VectorSIMD& operator -= (VectorSIMD<Real,12> a){
			_avx0=  _mm256_sub_pd(_avx0,a._avx0);
			_avx1=  _mm256_sub_pd(_avx1,a._avx1);
			_avx2=  _mm256_sub_pd(_avx2,a._avx2);
			return *this;
		}

		VectorSIMD operator -() const{
			const static __m256d m1 = _mm256_set_pd(-1.0,-1.0,-1.0,-1.0);
			VectorSIMD<double,12> result;
			result._avx0=  _mm256_mul_pd(this->_avx0,m1);
			result._avx1=  _mm256_mul_pd(this->_avx1,m1);
			result._avx2=  _mm256_mul_pd(this->_avx2,m1);
			return result;
		}

		void load(Real const* a){
			_avx0 =  _mm256_loadu_pd(&a[0]);
			_avx1 =  _mm256_loadu_pd(&a[4]);
			_avx2 =  _mm256_loadu_pd(&a[8]);
		}

		void load_aligned(Real const* a){
			_avx0 =  _mm256_load_pd(&a[0]);
			_avx1 =  _mm256_load_pd(&a[4]);
			_avx2 =  _mm256_load_pd(&a[8]);
		}

		void convert(Real *a) const {
			_mm256_storeu_pd(&a[0],_avx0);
			_mm256_storeu_pd(&a[4],_avx1);
			_mm256_storeu_pd(&a[8],_avx2);
		}

		void convert_aligned(Real *a) const {
			_mm256_storeu_pd(&a[0],_avx0);
			_mm256_storeu_pd(&a[4],_avx1);
			_mm256_storeu_pd(&a[8],_avx2);
		}

	};


	inline VectorSIMD<double,12> operator*(double a, VectorSIMD<double,12> b){
		VectorSIMD<double,12> c;
		__m256d _a =  _mm256_set_pd(a,a,a,a);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		c._avx1=  _mm256_mul_pd(_a, b._avx1);
		c._avx2=  _mm256_mul_pd(_a, b._avx2);
		return c;
	}


	inline VectorSIMD<double,12> operator*(VectorSIMD<double,12> a ,double b){
		VectorSIMD<double,12> c;
		__m256d _b =  _mm256_set_pd(b,b,b,b);
		c._avx0=  _mm256_mul_pd(a._avx0, _b);
		c._avx1=  _mm256_mul_pd(a._avx1, _b);
		c._avx2=  _mm256_mul_pd(a._avx2, _b);
		return c;
	}


	inline VectorSIMD<double,12> operator*(int a, VectorSIMD<double,12> b){
		if(a==1){return b;}
		VectorSIMD<double,12> c;
		double q = static_cast<double>(a);
		__m256d _a =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_a, b._avx0);
		c._avx1=  _mm256_mul_pd(_a, b._avx1);
		c._avx2=  _mm256_mul_pd(_a, b._avx2);
		return c;
	}


	inline VectorSIMD<double,12> operator*(VectorSIMD<double,12> a, int b){
		if(b==1){return a;}
		VectorSIMD<double,12> c;
		double q = static_cast<double>(b);
		__m256d _b =  _mm256_set_pd(q,q,q,q);
		c._avx0=  _mm256_mul_pd(_b, a._avx0);
		c._avx1=  _mm256_mul_pd(_b, a._avx1);
		c._avx2=  _mm256_mul_pd(_b, a._avx2);
		return c;
	}


	inline VectorSIMD<double,12> operator*(VectorSIMD<double,12> a, VectorSIMD<double,12> b){
		VectorSIMD<double,12> c;
		c._avx0=  _mm256_mul_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_mul_pd(a._avx1, b._avx1);
		c._avx2=  _mm256_mul_pd(a._avx2, b._avx2);
		return c;
	}


	inline VectorSIMD<double,12> operator+(VectorSIMD<double,12> a, VectorSIMD<double,12> b){
		VectorSIMD<double,12> c;
		c._avx0=  _mm256_add_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_add_pd(a._avx1, b._avx1);
		c._avx2=  _mm256_add_pd(a._avx2, b._avx2);
		return c;
	}


	inline VectorSIMD<double,12> operator-(VectorSIMD<double,12> a, VectorSIMD<double,12> b){
		VectorSIMD<double,12> c;
		c._avx0=  _mm256_sub_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_sub_pd(a._avx1, b._avx1);
		c._avx2=  _mm256_sub_pd(a._avx2, b._avx2);
		return c;
	}


	inline VectorSIMD<double,12> operator/(VectorSIMD<double,12> a, VectorSIMD<double,12> b){
		VectorSIMD<double,12> c;
		c._avx0=  _mm256_div_pd(a._avx0, b._avx0);
		c._avx1=  _mm256_div_pd(a._avx1, b._avx1);
		c._avx2=  _mm256_div_pd(a._avx2, b._avx2);
		return c;
	}


	inline VectorSIMD<double,12> exp(VectorSIMD<double,12> a){
		double a_d[12]; a.convert(a_d);
		a_d[0] = std::exp(a_d[0]);
		a_d[1] = std::exp(a_d[1]);
		a_d[2] = std::exp(a_d[2]);
		a_d[3] = std::exp(a_d[3]);
		a_d[4] = std::exp(a_d[4]);
		a_d[5] = std::exp(a_d[5]);
		a_d[6] = std::exp(a_d[6]);
		a_d[7] = std::exp(a_d[7]);
		a_d[8] = std::exp(a_d[8]);
		a_d[9] = std::exp(a_d[9]);
		a_d[10] = std::exp(a_d[10]);
		a_d[11] = std::exp(a_d[11]);
		VectorSIMD<double,12> r(a_d);
		return r;
	}

	inline VectorSIMD<double,12> sqrt(VectorSIMD<double,12> a){
		double a_d[12]; a.convert(a_d);
		a_d[0] = std::sqrt(a_d[0]);
		a_d[1] = std::sqrt(a_d[1]);
		a_d[2] = std::sqrt(a_d[2]);
		a_d[3] = std::sqrt(a_d[3]);
		a_d[4] = std::sqrt(a_d[4]);
		a_d[5] = std::sqrt(a_d[5]);
		a_d[6] = std::sqrt(a_d[6]);
		a_d[7] = std::sqrt(a_d[7]);
		a_d[8] = std::sqrt(a_d[8]);
		a_d[9] = std::sqrt(a_d[9]);
		a_d[10] = std::sqrt(a_d[10]);
		a_d[11] = std::sqrt(a_d[11]);
		VectorSIMD<double,12> r(a_d);
		return r;
	}

	inline std::ostream& operator<<(std::ostream& os, VectorSIMD<double,12> a){ 
		double ad[12];
		a.convert(ad);
		os << "{" << ad[0] << ","  << ad[1] << ","  << ad[2] << ","  << ad[3] << ","  << ad[4] << ","  << ad[5] << ","  << ad[6] << ","  << ad[7] << ","  << ad[8] << ","  << ad[9] << ","  << ad[10] << ","  << ad[11] << "}";
		return os;
	}  



	inline double reduce(VectorSIMD<double,12> a){
		__m256d s = _mm256_hadd_pd(a._avx0,a._avx0);
		double new_eri_double =  ((double*)&s)[0] + ((double*)&s)[2];
		s = _mm256_hadd_pd(a._avx1,a._avx1);
		new_eri_double +=  ((double*)&s)[0] + ((double*)&s)[2];
		s = _mm256_hadd_pd(a._avx2,a._avx2);
		new_eri_double +=  ((double*)&s)[0] + ((double*)&s)[2];
		return new_eri_double;
	}	
