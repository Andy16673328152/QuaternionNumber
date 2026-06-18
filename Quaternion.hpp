/*
Copyright 2026-----aleph_0 Andy16673328152(https://www.luogu.com.cn/user/1476867) All rights reserved.
Quaternion.hpp
This library provides quaternion basic operations,transcendental functions and literals.
*/
#ifndef QUATERNION_HPP
#define QUATERNION_HPP
#include <ios>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <tuple>
#include <limits>
#include <array>
#include <complex>
#if __cplusplus >= 202600L
#include <linalg>
#endif
#if __has_include(<glm/glm.hpp>)
#define HAS_GLM 1
#include <glm/glm.hpp>
#else
#define HAS_GLM 0
#endif
#if __has_include(<Eigen/Eigen>)
#define HAS_EIGEN 1
#include <Eigen/Eigen>
#else
#define HAS_EIGEN 0
#endif
#if __cplusplus >=202002L
#define MAYBE_CONSTEXPR constexpr
#include <numbers>
#else
#define MAYBE_CONSTEXPR inline
#endif
#ifdef __SSE__
#include <xmmintrin.h>
#endif
#ifdef __AVX2__
#include <immintrin.h>
#endif
//Constants
namespace Constants{
	#if __cplusplus>=202002L
	template<typename T>
	inline constexpr T PI=std::numbers::pi_v<T>;
	#else
	template<typename T>
	inline constexpr T PI=T(3.14159265358979323846);
	#endif
	template<typename T>
	inline constexpr T eps=T(1e-10);
}
template <typename T>
T stot(const std::string &s){
	static_assert(std::is_same<T,float>::value||std::is_same<T,double>::value||std::is_same<T,long double>::value,"T must be 'float','double' or 'long double',not 'int'");
	if constexpr(std::is_same<T,float>::value){
		return std::stof(s);
	}
	if constexpr(std::is_same<T,double>::value){
		return std::stod(s);
	}
	if constexpr(std::is_same<T,long double>::value){
		return std::stold(s);
	}
}
template <typename T>
struct alignas(16) Quaternion{
	static_assert(std::is_same<T,float>::value||std::is_same<T,double>::value||std::is_same<T,long double>::value,"T must be 'float','double' or 'long double',not 'int'");
	T r,i,j,k;
	constexpr T w()const noexcept{return r;}
	constexpr T x()const noexcept{return i;}
	constexpr T y()const noexcept{return j;}
	constexpr T z()const noexcept{return k;}
	constexpr Quaternion conj()const{
		return Quaternion(r,-i,-j,-k);
	}
	//getters and setters
	constexpr T real()const noexcept{
		return r;
	}
	void real(T val){
		r=val;
	}
	constexpr T imagi()const noexcept{
		return i;
	}
	void imagi(T val){
		i=val;
	}
	constexpr T imagj()const noexcept{
		return j;
	}
	void imagj(T val){
		j=val;
	}
	constexpr T imagk()const noexcept{
		return k;
	}
	void imagk(T val){
		k=val;
	}
	//constructors
	constexpr Quaternion()noexcept:r(0.0),i(0.0),j(0.0),k(0.0){}
	explicit constexpr Quaternion(T a,T b,T c,T d)noexcept:r(a),i(b),j(c),k(d){}
	explicit constexpr Quaternion(T a)noexcept:r(a){}
	explicit constexpr Quaternion(const std::complex<T> &a)noexcept:r(a.real()),i(a.imag()),j(0),k(0){}
	explicit constexpr Quaternion(const std::array<T,4> &a)noexcept:r(a[0]),i(a[1]),j(a[2]),k(a[3]){}
	explicit constexpr Quaternion(const std::array<T,3> &a)noexcept:r(0),i(a[0]),j(a[1]),k(a[2]){}
	explicit constexpr Quaternion(const std::initializer_list<T> &l)noexcept{
		if(l.size()==3){
			r=0;
			i=l.begin()[0];
			j=l.begin()[1];
			k=l.begin()[2];
		}else{
			if(l.size()==4){
				r=l.begin()[0];
				i=l.begin()[1];
				j=l.begin()[2];
				k=l.begin()[3];
			}else{
				r=i=j=k=0;
			}
		}
	}
	template <typename U>
	explicit Quaternion(const Quaternion<U> &a)noexcept:r(a.r),i(a.i),j(a.j),k(a.k){}
	Quaternion& operator =(const Quaternion &a){
		r=a.r;
		i=a.i;
		j=a.j;
		k=a.k;
		return *this;
	}
	template <typename U>
	Quaternion& operator =(const Quaternion<U> &a){
		r=a.r;
		i=a.i;
		j=a.j;
		k=a.k;
		return *this;
	}
	Quaternion& operator +=(T b){
		r+=b;
		return *this;
	}
	Quaternion& operator -=(T b){
		r-=b;
		return *this;
	}
	Quaternion& operator *=(T b){
		r*=b;
		i*=b;
		j*=b;
		k*=b;
		return *this;
	}
	Quaternion& operator /=(T b){
		r/=b;
		i/=b;
		j/=b;
		k/=b;
		return *this;
	}
	bool is_real()const{
		return (std::abs(i)<Constants::eps<T>)&&(std::abs(j)<Constants::eps<T>)&&(std::abs(k)<Constants::eps<T>);
	}
	bool is_imaginary()const{
		return (std::abs(r)<Constants::eps<T>);
	}
	template<std::size_t N>
	decltype(auto) get()const{
		static_assert(N<4,"Index out of range");
		if constexpr(N==0){
			return r;
		}
		if constexpr(N==1){
			return i;
		}
		if constexpr(N==2){
			return j;
		}
		if constexpr(N==3){
			return k;
		}
	}
	static Quaternion identity(){
		return Quaternion(T(1));
	}
	static Quaternion zero(){
		return Quaternion();
	}
	static Quaternion unit(T r,T i,T j,T k){
		Quaternion q(r,i,j,k);
		return q.normalized();
	}
	Quaternion& normalize(){
		T al=std::sqrt(r*r+i*i+j*j+k*k);
		if(al>Constants::eps<T>){
			T ial=T(1)/al;
			if constexpr(std::is_same<T,float>::value){
				#ifdef __SSE__
				__m128 vq=_mm_load_ps(&(this->r));
				__m128 fourxial=_mm_set_ps(ial,ial,ial,ial);
				__m128 rr=_mm_mul_ps(vq,fourxial);
				float rt[4];
				_mm_store_ps(rt,rr);
				r=rt[0];
				i=rt[1];
				j=rt[2];
				k=rt[3];
				#else
				r*=ial;
				i*=ial;
				j*=ial;
				k*=ial;
				#endif
			}else{
				if constexpr(std::is_same<T,double>::value){
					#ifdef __AVX2__
					__m256d vq=_mm256_load_pd(&(this->r));
					__m256d fourxial=_mm256_set_pd(ial,ial,ial,ial);
					__m256d rr=_mm256_mul_pd(vq,fourxial);
					double rt[4];
					_mm256_store_pd(rt,rr);
					r=rt[0];
					i=rt[1];
					j=rt[2];
					k=rt[3];
					#else
					r*=ial;
					i*=ial;
					j*=ial;
					k*=ial;
					#endif
				}else{
					r*=ial;
					i*=ial;
					j*=ial;
					k*=ial;	
				}
			}
		}
		return *this;
	}
	Quaternion normalized()const{
		T al=std::sqrt(r*r+i*i+j*j+k*k);
		if(al>Constants::eps<T>){
			T ial=T(1)/al;
			if constexpr(std::is_same<T,float>::value){
				#ifdef __SSE__
				__m128 vq=_mm_load_ps(&(this->r));
				__m128 fourxial=_mm_set_ps(ial,ial,ial,ial);
				__m128 rr=_mm_mul_ps(vq,fourxial);
				float rt[4];
				_mm_store_ps(rt,rr);
				return Quaternion(rt[0],rt[1],rt[2],rt[3]);
				#else
				return Quaternion(r*ial,i*ial,j*ial,k*ial);
				#endif
			}else{
				if constexpr(std::is_same<T,double>::value){
					#ifdef __AVX2__
					__m256d vq=_mm256_load_pd(&(this->r));
					__m256d fourxial=_mm256_set_pd(ial,ial,ial,ial);
					__m256d rr=_mm256_mul_pd(vq,fourxial);
					double rt[4];
					_mm256_store_pd(rt,rr);
					return Quaternion(rt[0],rt[1],rt[2],rt[3]);
					#else
					return Quaternion(r*ial,i*ial,j*ial,k*ial);
					#endif
				}else{
					return Quaternion(r*ial,i*ial,j*ial,k*ial);	
				}
			}
		}
		return *this;
	}
	Quaternion normalized(T norm)const{
		T al=std::sqrt(norm);
		if(al>Constants::eps<T>){
			T ial=T(1)/al;
			if constexpr(std::is_same<T,float>::value){
				#ifdef __SSE__
				__m128 vq=_mm_load_ps(&(this->r));
				__m128 fourxial=_mm_set_ps(ial,ial,ial,ial);
				__m128 r=_mm_mul_ps(vq,fourxial);
				float rt[4];
				_mm_store_ps(rt,r);
				return Quaternion(rt[0],rt[1],rt[2],rt[3]);
				#else
				return Quaternion(r*ial,i*ial,j*ial,k*ial);
				#endif
			}else{
				if constexpr(std::is_same<T,double>::value){
					#ifdef __AVX2__
					__m256d vq=_mm256_load_pd(&(this->r));
					__m256d fourxial=_mm256_set_pd(ial,ial,ial,ial);
					__m256d r=_mm256_mul_pd(vq,fourxial);
					double rt[4];
					_mm256_store_pd(rt,r);
					return Quaternion(rt[0],rt[1],rt[2],rt[3]);
					#else
					return Quaternion(r*ial,i*ial,j*ial,k*ial);
					#endif
				}else{
					return Quaternion(r*ial,i*ial,j*ial,k*ial);	
				}
			}
		}
		return *this;
	}
	//Determine whether a quaternion number is a unit quaternion number
	MAYBE_CONSTEXPR bool isunit(T eeps=Constants::eps<T>)const{
		T alsq=r*r+i*i+j*j+k*k;
		return std::abs(alsq-1)<eeps;
	}
	//(x,y,z),theta->quaternion
	MAYBE_CONSTEXPR static Quaternion FromAxisAngle(T angle,T x,T y,T z){
		T semiangle=angle/2.0;
		Quaternion q(std::cos(semiangle),x*std::sin(semiangle),y*std::sin(semiangle),z*std::sin(semiangle));
		return q.normalized();
	} 
	//quaternion->rotation matrix
	MAYBE_CONSTEXPR std::array<T,9> ToRotationMatrix()const{
		T ww=r*r,xx=i*i,yy=j*j,zz=k*k,xw=r*i,yw=r*j,zw=r*k,xy=i*j,xz=i*k,yz=j*k;
		return {
			1-2*(yy+zz),2*(xy-zw),2*(xz+yw),
			2*(xy+zw),1-2*(xx+zz),2*(yz-xw),
			2*(xz-yw),2*(yz+xw),1-2*(xx+yy)
		};
	}
	//rotation matrix->quaternion
	static Quaternion FromRotationMatrix(const std::array<T,9>& mat){
		T t=mat[0]+mat[4]+mat[8];
		Quaternion q;
		if(t>0){
	 		T s=2.0*(std::sqrt(t+1));
	 		T is=T(1)/s;
	 		q.r=0.25*s;
	 		q.i=(mat[5]-mat[7])*is;
	 		q.j=(mat[6]-mat[2])*is;
	 		q.k=(mat[1]-mat[3])*is;
		}else{
		 	if(mat[0]>mat[4]&&mat[0]>mat[8]){
		 		T s=2.0*std::sqrt(T(1)+mat[0]-mat[4]-mat[8]);
		 		T is=T(1)/s;
		 		q.r=(mat[5]-mat[7])*is;
		 		q.i=0.25*s;
		 		q.j=(mat[3]+mat[1])*is;
		 		q.k=(mat[6]+mat[2])*is;
			}else{
				if(mat[4]>mat[8]){
					T s=2.0*std::sqrt(T(1)+mat[4]-mat[0]-mat[8]);
					q.r=(mat[6]-mat[2])/s;
					q.i=(mat[3]+mat[1])/s;
					q.j=0.25*s;
					q.k=(mat[7]+mat[5])/s;
				}else{
					T s=2.0*std::sqrt(T(1)+mat[8]-mat[4]-mat[0]);
					T is=T(1)/s;
					q.r=(mat[1]-mat[3])*is;
					q.i=(mat[6]+mat[2])*is;
					q.j=(mat[7]+mat[5])*is;
					q.k=0.25*s;
				}
			}
		}
		return q;
	}
	enum class EulerOrder{
		XYZ,XZY,YXZ,YZX,ZXY,ZYX
	};
	//quaternion->(roll,pitch,yaw)
	MAYBE_CONSTEXPR std::array<T,3> ToEulerAngles(EulerOrder order=EulerOrder::XYZ,T epss=T(1)-Constants::eps<T>)const{
		Quaternion q=this->normalized();
		T w=q.r,x=q.i,y=q.j,z=q.k;
		T roll,pitch,yaw;
		auto clasin=[](T val){
			return std::asin(std::max(T(-1),std::min(T(1),val)));
		};
		switch(order){
			case EulerOrder::XYZ:
				pitch=clasin(2*(w*y-z*x));
				if(std::abs(std::sin(pitch))<epss){
					roll=std::atan2(2*(w*x+y*z),1-2*(x*x+y*y));
					yaw=std::atan2(2*(w*z+x*y),1-2*(y*y+z*z));
				}else{
					roll=std::atan2(2*(x*z-w*y),1-2*(x*x+z*z));
					yaw=0;
				}
				return {roll,pitch,yaw};
			case EulerOrder::XZY:
				yaw=clasin(2*(w*z+x*y));
				if(std::abs(std::sin(yaw))<epss){
					pitch=std::atan2(2*(w*y-x*z),1-2*(y*y+z*z));
					roll=std::atan2(2*(w*x-y*z),1-2*(x*x+z*z));
				}else{
					roll=std::atan2(2*(x*z+w*y),1-2*(x*x+y*y));
					pitch=0;
				}
				return {roll,pitch,yaw};
			case EulerOrder::YXZ:
				roll=clasin(2*(w*x-y*z));
				if(std::abs(std::sin(roll))<epss){
					pitch=std::atan2(2*(w*y+x*z),1-2*(x*x+y*y));
					yaw=std::atan2(2*(w*z+x*y),1-2*(x*x+z*z));
				}else{
					yaw=std::atan2(2*(x*z+w*y),1-2*(y*y+z*z));
					pitch=0;
				}
				return {roll,pitch,yaw};
			case EulerOrder::YZX:
				yaw=clasin(2*(w*z-x*y));
				if(std::abs(std::sin(yaw))<epss){
					pitch=std::atan2(2*(w*x+y*z),1-2*(x*x+z*z));
					roll=std::atan2(2*(w*y+x*z),1-2*(y*y+z*z));
				}else{
					roll=std::atan2(2*(x*y+w*z),1-2*(x*x+y*y));
					pitch=0;
				}
				
				return {roll,pitch,yaw};
			case EulerOrder::ZXY:
				roll=clasin(2*(w*x+y*z));
				if(std::abs(std::sin(roll))<epss){
					pitch=std::atan2(2*(w*y-x*z),1-2*(x*x+y*y));
					yaw=std::atan2(2*(w*z-x*y),1-2*(x*x+z*z));
				}else{
					yaw=std::atan2(2*(x*y+w*z),1-2*(y*y+z*z));
					pitch=0;
				}
				
				return {roll,pitch,yaw};
			case EulerOrder::ZYX:
				pitch=clasin(2*(w*y-z*x));
				if(std::abs(std::sin(pitch))<epss){
					roll=std::atan2(2*(w*x+y*z),1-2*(x*x+y*y));
					yaw=std::atan2(2*(w*z+x*y),1-2*(y*y+z*z));
				}else{
					yaw=std::atan2(2*(x*z-w*y),1-2*(x*x+z*z));
					roll=0;
				}
				return {roll,pitch,yaw};
			default:
				return {0,0,0};
		}
	}
	//(roll,pitch,yaw)->quaternion
	static Quaternion FromEulerAngles(T roll,T pitch,T yaw,EulerOrder order=EulerOrder::XYZ){
		T c1=std::cos(roll/2);
		T s1=std::sin(roll/2);
		T c2=std::cos(pitch/2);
		T s2=std::sin(pitch/2);
		T c3=std::cos(yaw/2);
		T s3=std::sin(yaw/2);
		Quaternion q;
		switch(order){
			case EulerOrder::XYZ:
				q.r=c1*c2*c3-s1*s2*s3;
				q.i=s1*c2*c3+c1*s2*s3;
				q.j=c1*s2*c3-s1*c2*s3;
				q.k=c1*c2*s3+s1*s2*c3;
				return q.normalized();
			case EulerOrder::XZY:
				q.r=c1*c2*c3+s1*s2*s3;
				q.i=s1*c2*c3-c1*s2*s3;
				q.j=c1*c2*s3-s1*s2*c3;
				q.k=c1*s2*c3+s1*c2*s3;
				return q.normalized();
			case EulerOrder::YXZ:
				q.r=c1*c2*c3+s1*s2*s3;
				q.i=s1*c2*c3+c1*s2*s3;
				q.j=c1*s2*c3-s1*c2*s3;
				q.k=c1*c2*s3-s1*s2*c3;
				return q.normalized();
			case EulerOrder::YZX:
				q.r=c1*c2*c3-s1*s2*s3;
				q.i=c1*s2*c3+s1*c2*s3;
				q.j=s1*c2*c3+c1*s2*s3;
				q.k=c1*c2*s3-s1*s2*c3;
				return q.normalized();
			case EulerOrder::ZXY:
				q.r=c1*c2*c3-s1*s2*s3;
				q.i=c1*s2*c3-s1*c2*s3;
				q.j=c1*c2*s3+s1*s2*c3;
				q.k=s1*c2*c3+c1*s2*s3;
				return q.normalized();
			case EulerOrder::ZYX:
				q.r=c1*c2*c3+s1*s2*s3;
				q.i=c1*c2*s3-s1*s2*c3;
				q.j=c1*s2*c3+s1*c2*s3;
				q.k=s1*c2*c3-c1*s2*s3;
				return q.normalized();
			default:
				return Quaternion(1,0,0,0);
		}
	}
	MAYBE_CONSTEXPR std::array<T,3> RotateVector(T x,T y,T z)const{
		T rx=2*(j*z-k*y),ry=2*(k*x-i*z),rz=2*(i*y-j*x);
		return {x+r*rx+(j*rz-k*ry),y+r*ry+(k*rx-i*rz),z+r*rz+(i*ry-j*rx)};
	}
	MAYBE_CONSTEXPR std::array<T,3> RotateVector(const std::array<T,3> &v)const{
		return RotateVector(v[0],v[1],v[2]);
	}
	//axis(x,y,z)
	MAYBE_CONSTEXPR std::array<T,3> axis()const{
		if(std::isnan(r)||std::isnan(i)||std::isnan(j)||std::isnan(k)){
			return {1,0,0};
		}
		T cr=std::max(T(-1),std::min(T(1),r));
		T a=std::sqrt(T(1)-cr*cr);
		if(a<Constants::eps<T>||std::isnan(a)){
			return {1,0,0};
		}
		T ia=T(1)/a;
		return {i*ia,j*ia,k*ia};
	}
	//theta
	MAYBE_CONSTEXPR T angle()const{
		if(std::isnan(r)||std::isnan(i)||std::isnan(j)||std::isnan(k)){
			return T(0);
		}
		T cr=std::max(T(-1),std::min(T(1),r));
		return T(2)*std::acos(cr);
	}
	MAYBE_CONSTEXPR T angleDegrees()const{
		return angle()*T(180)/Constants::PI<T>;
	}
	//T(1)/(*this)
	MAYBE_CONSTEXPR Quaternion inverse()const{
		if(iszero(*this)){
			return Quaternion();
		}
		return (this->conj())/(r*r+i*i+j*j+k*k);
	}
};
//tuple
namespace std{
	template<typename T>
	struct tuple_size<Quaternion<T> >:std::integral_constant<size_t,4>{};
	template<typename T>
	struct tuple_element<0,Quaternion<T> >{using type=T;};
	template<typename T>
	struct tuple_element<1,Quaternion<T> >{using type=T;};
	template<typename T>
	struct tuple_element<2,Quaternion<T> >{using type=T;};
	template<typename T>
	struct tuple_element<3,Quaternion<T> >{using type=T;};
}
#if __cplusplus >=202600L
template<typename T>
std::linalg::quaternion<T> to_std(const Quaternion<T> &q){
	return std::linalg::quaternion<T>(q.r,q.i,q.j,q.k);
}
template<typename T>
Quaternion<T> from_std(const std::linalg::quaternion<T> &q){
	return Quaternion<T>(q.w(),q.x(),q.y(),q.z());
}
#endif
#if HAS_GLM==1
glm::quat to_glmf(const Quaternion<float> &q){
	glm::quat g;
	g.w=q.r;
	g.x=q.i;
	g.y=q.j;
	g.z=q.k;
	return g;
}
glm::dquat to_glmd(const Quaternion<double> &q){
	glm::dquat g;
	g.w=q.r;
	g.x=q.i;
	g.y=q.j;
	g.z=q.k;
	return g;
}
Quaternion<float> from_glmf(const glm::quat &g){
	return Quaternion<float>(g.w,g.x,g.y,g.z);
}
Quaternion<double> from_glmd(const glm::dquat &g){
	return Quaternion<double>(g.w,g.x,g.y,g.z);
}
#endif
#if HAS_EIGEN==1
Quaternionf to_eigenqf(const Quaternion<float> &q){
	return Quaternionf(q.r,q.i,q.j,q.k);
}
Quaterniond to_eigenqd(const Quaternion<double> &q){
	return Quaterniond(q.r,q.i,q.j,q.k);
}
Quaternion<float> from_eigenqf(const Quaternionf &q){
	return Quaternion<float>(q.w(),q.x(),q.y(),q.z());
}
Quaternion<double> from_eigenqd(const Quaterniond &q){
	return Quaternion<double>(q.w(),q.x(),q.y(),q.z());
}
#endif
//basic operators
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> operator +(const Quaternion<T> &a,T b){
	return Quaternion<T>(a.r+b,a.i,a.j,a.k);
} 
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> operator +(T a,const Quaternion<T> &b){
	return Quaternion<T>(b.r+a,b.i,b.j,b.k);
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> operator -(const Quaternion<T> &a,T b){
	return Quaternion<T>(a.r-b,a.i,a.j,a.k);
} 
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> operator -(const Quaternion<T> &a){
	return Quaternion<T>(-a.r,-a.i,-a.j,-a.k);
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> operator -(T a,const Quaternion<T> &b){
	return Quaternion<T>(a-b.r,-b.i,-b.j,-b.k);
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> operator *(const Quaternion<T> &a,T b){
	return Quaternion<T>(a.r*b,a.i*b,a.j*b,a.k*b);
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> operator *(T a,const Quaternion<T> &b){
	return Quaternion<T>(b.r*a,b.i*a,b.j*a,b.k*a);
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> operator /(Quaternion<T> a,T b){
	if(std::abs(b)<Constants::eps<T>){
		throw std::runtime_error("Division by zero");
	}
	T ib=T(1)/b;
	return Quaternion<T>(a.r*ib,a.i*ib,a.j*ib,a.k*ib);
}
template<typename T>
MAYBE_CONSTEXPR T abs(const Quaternion<T> &a){
	return std::sqrt(a.r*a.r+a.i*a.i+a.j*a.j+a.k*a.k);
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> operator /(T a,const Quaternion<T> &b){
	T mb=abs(b);
	if(mb<Constants::eps<T>){
		throw std::runtime_error("Division by zero");
	}
	return a*(b.conj()/(mb*mb));
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> operator +(const Quaternion<T> &a,const Quaternion<T> &b){
	if constexpr(std::is_same<T,float>::value){
	#ifdef __SSE__
	__m128 va=_mm_load_ps(&a.r);
	__m128 vb=_mm_load_ps(&b.r);
	__m128 r=_mm_add_ps(va,vb);
	float rs[4];
	_mm_store_ps(rs,r);
	return Quaternion<T>(rs[0],rs[1],rs[2],rs[3]);
	#else
	return Quaternion<T>(a.r+b.r,a.i+b.i,a.j+b.j,a.k+b.k);
	#endif
	}else{
		if constexpr(std::is_same<T,double>::value){
			#ifdef __AVX2__
			__m256d va=_mm256_load_pd(&a.r);
			__m256d vb=_mm256_load_pd(&b.r);
			__m256d r=_mm256_add_pd(va,vb);
			double rs[4];
			_mm256_store_pd(rs,r);
			return Quaternion<T>(rs[0],rs[1],rs[2],rs[3]);
			#else
			return Quaternion<T>(a.r+b.r,a.i+b.i,a.j+b.j,a.k+b.k);
			#endif
		}else{
			return Quaternion<T>(a.r+b.r,a.i+b.i,a.j+b.j,a.k+b.k);
		}
	}
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> operator -(const Quaternion<T> &a,const Quaternion<T> &b){
	if constexpr(std::is_same<T,float>::value){
	#ifdef __SSE__
	__m128 va=_mm_load_ps(&a.r);
	__m128 vb=_mm_load_ps(&b.r);
	__m128 r=_mm_sub_ps(va,vb);
	float rs[4];
	_mm_store_ps(rs,r);
	return Quaternion<T>(rs[0],rs[1],rs[2],rs[3]);
	#else
	return Quaternion<T>(a.r-b.r,a.i-b.i,a.j-b.j,a.k-b.k);
	#endif
	}else{
		if constexpr(std::is_same<T,double>::value){
			#ifdef __AVX2__
			__m256d va=_mm256_load_pd(&a.r);
			__m256d vb=_mm256_load_pd(&b.r);
			__m256d r=_mm256_sub_pd(va,vb);
			double rs[4];
			_mm256_store_pd(rs,r);
			return Quaternion<T>(rs[0],rs[1],rs[2],rs[3]);
			#else
			return Quaternion<T>(a.r-b.r,a.i-b.i,a.j-b.j,a.k-b.k);
			#endif
		}else{
			return Quaternion<T>(a.r-b.r,a.i-b.i,a.j-b.j,a.k-b.k);
		}
	}
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> operator *(Quaternion<T> a,Quaternion<T> b){
	if constexpr(std::is_same<T,float>::value){
		#ifdef __SSE__
		__m128 va=_mm_load_ps(&a.r);
		__m128 vb=_mm_load_ps(&b.r);
		float tp[4];
		_mm_store_ps(tp,vb);
		__m128 a0=_mm_shuffle_ps(va,va,_MM_SHUFFLE(0,0,0,0));
		__m128 a1=_mm_shuffle_ps(va,va,_MM_SHUFFLE(1,1,1,1));
		__m128 a2=_mm_shuffle_ps(va,va,_MM_SHUFFLE(2,2,2,2));
		__m128 a3=_mm_shuffle_ps(va,va,_MM_SHUFFLE(3,3,3,3));
		__m128 r=_mm_mul_ps(a0,_mm_set_ps(tp[3],tp[2],tp[1],tp[0]));
		r=_mm_add_ps(r,_mm_mul_ps(a1,_mm_set_ps(tp[2],-tp[3],tp[0],-tp[1])));
		r=_mm_add_ps(r,_mm_mul_ps(a2,_mm_set_ps(-tp[1],tp[0],tp[3],-tp[2])));
		r=_mm_add_ps(r,_mm_mul_ps(a3,_mm_set_ps(tp[0],tp[1],-tp[2],-tp[3])));
		float rs[4];
		_mm_store_ps(rs,r);
		return Quaternion<T>(rs[0],rs[1],rs[2],rs[3]);
		#else
		return Quaternion<T>(a.r*b.r-a.i*b.i-a.j*b.j-a.k*b.k,a.r*b.i+a.i*b.r+a.j*b.k-a.k*b.j,a.r*b.j+a.j*b.r+a.k*b.i-a.i*b.k,a.r*b.k+a.k*b.r+a.i*b.j-a.j*b.i);
		#endif
	}else{
		if constexpr(std::is_same<T,double>::value){
			#ifdef __AVX2__
			__m256d va=_mm256_load_pd(&a.r);
			__m256d vb=_mm256_load_pd(&b.r);
			double tp[4];
			_mm256_store_pd(tp,vb);
			__m256d a0=_mm256_permute4x64_pd(va,_MM_SHUFFLE(0,0,0,0));
			__m256d a1=_mm256_permute4x64_pd(va,_MM_SHUFFLE(1,1,1,1));
			__m256d a2=_mm256_permute4x64_pd(va,_MM_SHUFFLE(2,2,2,2));
			__m256d a3=_mm256_permute4x64_pd(va,_MM_SHUFFLE(3,3,3,3));
			__m256d r=_mm256_mul_pd(a0,vb);
			r=_mm256_add_pd(r,_mm256_mul_pd(a1,_mm256_set_pd(tp[2],-tp[3],tp[0],-tp[1])));
			r=_mm256_add_pd(r,_mm256_mul_pd(a2,_mm256_set_pd(-tp[1],tp[0],tp[3],-tp[2])));
			r=_mm256_add_pd(r,_mm256_mul_pd(a3,_mm256_set_pd(tp[0],tp[1],-tp[2],-tp[3])));
			double rs[4];
			_mm256_store_pd(rs,r);
			return Quaternion<T>(rs[0],rs[1],rs[2],rs[3]);
			#else
			return Quaternion<T>(a.r*b.r-a.i*b.i-a.j*b.j-a.k*b.k,a.r*b.i+a.i*b.r+a.j*b.k-a.k*b.j,a.r*b.j+a.j*b.r+a.k*b.i-a.i*b.k,a.r*b.k+a.k*b.r+a.i*b.j-a.j*b.i);
			#endif
		}else{
			return Quaternion<T>(a.r*b.r-a.i*b.i-a.j*b.j-a.k*b.k,a.r*b.i+a.i*b.r+a.j*b.k-a.k*b.j,a.r*b.j+a.j*b.r+a.k*b.i-a.i*b.k,a.r*b.k+a.k*b.r+a.i*b.j-a.j*b.i);
		}
	}
}
template<typename T>
MAYBE_CONSTEXPR std::array<T,3> operator *(const Quaternion<T> &a,std::array<T,3> v){
	Quaternion<T> h=pure(a*Quaternion(v)*a.inverse());
	return {h.i,h.j,h.k};
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> operator /(const Quaternion<T> &a,const Quaternion<T> &b){
	return a*(T(1)/b);
}
template<typename T>
MAYBE_CONSTEXPR bool isnan(const Quaternion<T> &a){
	return std::isnan(a.r)||std::isnan(a.i)||std::isnan(a.j)||std::isnan(a.k);
}
template<typename T>
MAYBE_CONSTEXPR bool isinf(const Quaternion<T> &a){
	return std::isinf(a.r)||std::isinf(a.i)||std::isinf(a.j)||std::isinf(a.k);
}
template<typename T>
MAYBE_CONSTEXPR bool ispure(const Quaternion<T> &a){
	return std::abs(a.r)<Constants::eps<T>;
}
template<typename T>
MAYBE_CONSTEXPR bool isnormalized(const Quaternion<T> &a){
	return std::abs(abs(a)-T(1))<Constants::eps<T>; 
}
template<typename T>
MAYBE_CONSTEXPR bool iszero(const Quaternion<T> &a){
	return std::abs(a.r)<Constants::eps<T>&&std::abs(a.i)<Constants::eps<T>&&std::abs(a.j)<Constants::eps<T>&&std::abs(a.k)<Constants::eps<T>;
}
template<typename T>
MAYBE_CONSTEXPR bool isidentity(const Quaternion<T> &a){
	return std::abs(a.r-T(1))<Constants::eps<T>&&std::abs(a.i)<Constants::eps<T>&&std::abs(a.j)<Constants::eps<T>&&std::abs(a.k)<Constants::eps<T>;
}
template<typename T>
MAYBE_CONSTEXPR bool operator ==(const Quaternion<T> &a,const Quaternion<T> &b){
	if(isnan(a)||isnan(b)){
		return false;
	}
	if(isinf(a)||isinf(b)){
		return (a.r==b.r)&&(a.i==b.i)&&(a.j==b.j)&&(a.k==b.k);
	}
	return (std::abs(a.r-b.r)<Constants::eps<T>)&&(std::abs(a.i-b.i)<Constants::eps<T>)&&(std::abs(a.j-b.j)<Constants::eps<T>)&&(std::abs(a.k-b.k)<Constants::eps<T>);
}
template<typename T>
MAYBE_CONSTEXPR bool almost_equal(const Quaternion<T> &a,const Quaternion<T> &b,T eeps=Constants::eps<T>){
	if(isnan(a)||isnan(b)){
		return false;
	}
	if(isinf(a)||isinf(b)){
		return (a.r==b.r)&&(a.i==b.i)&&(a.j==b.j)&&(a.k==b.k);
	}
	return (std::abs(a.r-b.r)<eeps)&&(std::abs(a.i-b.i)<eeps)&&(std::abs(a.j-b.j)<eeps)&&(std::abs(a.k-b.k)<eeps);
}
template<typename T>
MAYBE_CONSTEXPR bool operator !=(const Quaternion<T> &a,const Quaternion<T> &b){
	return !(a==b);
}
template<typename T>
Quaternion<T>& operator +=(Quaternion<T> &a,const Quaternion<T> &b){
	a=a+b;
	return a;
}
template<typename T>
Quaternion<T>& operator -=(Quaternion<T> &a,const Quaternion<T> &b){
	a=a-b;
	return a;
}
template<typename T>
Quaternion<T>& operator *=(Quaternion<T> &a,const Quaternion<T> &b){
	a=a*b;
	return a;
}
template<typename T>
Quaternion<T>& operator /=(Quaternion<T> &a,const Quaternion<T> &b){
	a=a/b;
	return a;
}
//pure quaternion
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> pure(const Quaternion<T> &a){
	return Quaternion<T>(0.0,a.i,a.j,a.k);
}
//logarithmic functions
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> log(const Quaternion<T> &a){
	if(isnan(a)){
		return Quaternion<T>();
	}
	T mg=abs(a);
	T pa=std::sqrt(a.i*a.i+a.j*a.j+a.k*a.k);
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::log(a.r),0,0,0);
	}
	T ipa=T(1)/pa;
	T img=T(1)/mg;
	Quaternion<T> u(0.0,a.i*ipa,a.j*ipa,a.k*ipa);
	return std::log(mg)+u*std::acos(a.r*img);
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> log10(const Quaternion<T> &a){
	return log(a)/std::log(10);
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> log2(const Quaternion<T> &a){
	return log(a)/std::log(2);
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> logp(const Quaternion<T> &a,const Quaternion<T> &b){
	return log(a)/log(b);
}
//power functions
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> exp(const Quaternion<T> &a){
	if(isnan(a)){
		return Quaternion<T>();
	}
	T pa=std::sqrt(a.i*a.i+a.j*a.j+a.k*a.k);
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::exp(a.r),0,0,0);
	}
	T ipa=T(1)/pa;
	Quaternion<T> u(0.0,a.i*ipa,a.j*ipa,a.k*ipa);
	return std::exp(a.r)*(std::cos(pa)+u*std::sin(pa));
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> pow(const Quaternion<T> &a,const Quaternion<T> &b){
	return exp(log(a)*b);
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> nthrt(const Quaternion<T> &a,const Quaternion<T> &b){
	return pow(a,T(1)/b);
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> sqrt(const Quaternion<T> &a){
	return nthrt(a,Quaternion<T>(2.0));
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> cbrt(const Quaternion<T> &a){
	return nthrt(a,Quaternion<T>(3.0));
}
//a ·b 
template<typename T>
MAYBE_CONSTEXPR T dot(const Quaternion<T> &a,const Quaternion<T> &b){
	return a.r*b.r+a.i*b.i+a.j*b.j+a.k*b.k;
}

template <typename T>
Quaternion<T> unit(T b,T c,T d){
	return Quaternion<T>(0.0,b/std::sqrt(b*b+c*c+d*d),c/std::sqrt(b*b+c*c+d*d),d/std::sqrt(b*b+c*c+d*d));
}
//sphere linear interpolation
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> slerp(const Quaternion<T> &a,const Quaternion<T> &b,T t,T eeps=Constants::eps<T>){
	T d=dot(a,b),sign=T(1);
	if(d<0){
		d=-d;
		sign=-T(1);
	}
	const T dt=T(1)-eeps;
	if(d>dt){
		Quaternion<T> q=a+(b*sign-a)*t;
		return q.normalized();
	}
	T theta=std::acos(d)*t;
	T sin_theta=std::sin(theta);
	T sin_theta0=std::sqrt(1-d*d);
	T s2=sin_theta/sin_theta0;
	T s1=std::cos(theta)-d*s2;
	return (a*s1)+(b*sign*s2);
}
//n-th root
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> omega(const Quaternion<T> &a,int n,int k){
	if(isnan(a)){
		return Quaternion<T>();
	}
	k%=n;
	if(a.is_real()){
		if(std::abs(a.r)<Constants::eps<T>){
			return Quaternion(0.0);
		}else{
			if(a.r>0){
				T poww=std::pow(a.r,T(1)/n);
				Quaternion<T> q(std::cos(2*Constants::PI<T>*k/n),std::sin(2*Constants::PI<T>*k/n),0,0);
				return poww*q;
			}else{
				T poww=std::pow(-a.r,T(1)/n);
				Quaternion<T> q(std::cos((Constants::PI<T>+2*Constants::PI<T>*k)/n),std::sin((Constants::PI<T>+2*Constants::PI<T>*k)/n),0,0);
				return poww*q;
			}
		}
	}
	T r=abs(a),theta=a.angle();
	Quaternion<T> u=pure(a)/std::sqrt(a.i*a.i+a.j*a.j+a.k*a.k);
	return std::pow(r,T(1)/n)*(std::cos((theta+2*Constants::PI<T>*k)/n)+u*std::sin((theta+2*Constants::PI<T>*k)/n));
}
//trigonometric functions
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> sin(const Quaternion<T> &a){
	if(isnan(a)){
		return Quaternion<T>();
	}
	T pa=std::sqrt(a.i*a.i+a.j*a.j+a.k*a.k);
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::sin(a.r),0,0,0);
	}
	T ipa=T(1)/pa;
	Quaternion<T> u(0.0,a.i*ipa,a.j*ipa,a.k*ipa);
	return std::sin(a.r)*std::cosh(pa)+u*std::cos(a.r)*std::sinh(pa);
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> cos(const Quaternion<T> &a){
	if(isnan(a)){
		return Quaternion<T>();
	}
	T pa=std::sqrt(a.i*a.i+a.j*a.j+a.k*a.k);
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::cos(a.r),0,0,0);
	}
	T ipa=T(1)/pa;
	Quaternion<T> u(0.0,a.i*ipa,a.j*ipa,a.k*ipa);
	return std::cos(a.r)*std::cosh(pa)-u*std::sin(a.r)*std::sinh(pa);
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> tan(const Quaternion<T> &a){
	return sin(a)/cos(a);
}
//inverse trigonometric functions
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> asin(const Quaternion<T> &a){
	if(isnan(a)){
		return Quaternion<T>();
	}
	T pa=std::sqrt(a.i*a.i+a.j*a.j+a.k*a.k);
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::asin(a.r),0,0,0);
	}
	T ipa=T(1)/pa;
	Quaternion<T> u(0.0,a.i*ipa,a.j*ipa,a.k*ipa);
	return -u*log(a*u+sqrt(T(1)-a*a));
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> acos(const Quaternion<T> &a){
	if(isnan(a)){
		return Quaternion<T>();
	}
	T pa=std::sqrt(a.i*a.i+a.j*a.j+a.k*a.k);
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::acos(a.r),0,0,0);
	}
	T ipa=T(1)/pa;
	Quaternion<T> u(0.0,a.i*ipa,a.j*ipa,a.k*ipa);
	return -u*log(a+u*sqrt(T(1)-a*a));
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> atan(const Quaternion<T> &a){
	if(isnan(a)){
		return Quaternion<T>();
	}
	T pa=std::sqrt(a.i*a.i+a.j*a.j+a.k*a.k);
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::atan(a.r),0,0,0);
	}
	T ipa=T(1)/pa;
	Quaternion<T> u(0.0,a.i*ipa,a.j*ipa,a.k*ipa);
	return (u/2.0)*log((u+a)/(u-a));
}
//hyperbolic functions
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> sinh(const Quaternion<T> &a){
	if(isnan(a)){
		return Quaternion<T>();
	}
	T pa=std::sqrt(a.i*a.i+a.j*a.j+a.k*a.k);
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::sinh(a.r),0,0,0);
	}
	T ipa=T(1)/pa;
	Quaternion<T> u(0.0,a.i*ipa,a.j*ipa,a.k*ipa);
	return std::sinh(a.r)*std::cos(pa)+u*std::cosh(a.r)*std::sin(pa);
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> cosh(const Quaternion<T> &a){
	if(isnan(a)){
		return Quaternion<T>();
	}
	T pa=std::sqrt(a.i*a.i+a.j*a.j+a.k*a.k);
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::cosh(a.r),0,0,0);
	}
	T ipa=T(1)/pa;
	Quaternion<T> u(0.0,a.i*ipa,a.j*ipa,a.k*ipa);
	return std::cosh(a.r)*std::cos(pa)+u*std::sinh(a.r)*std::sin(pa);
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> tanh(const Quaternion<T> &a){
	return sinh(a)/cosh(a);
}
//inverse hyperbolic functions
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> asinh(const Quaternion<T> &a){
	if(isnan(a)){
		return Quaternion<T>();
	}
	return log(a+sqrt(a*a+T(1)));
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> acosh(const Quaternion<T> &a){
	if(isnan(a)){
		return Quaternion<T>();
	}
	if(isidentity(-a)){
		throw std::domain_error("acosh:singular point -1");
	}
	T w2=a.r*a.r,v2=abs(pure(a));
	if(std::abs(w2-v2-T(1))<Constants::eps<T>&&v2>Constants::eps<T>){
		Quaternion<T> sq=sqrt(T(1)-a*a);
		return log(a+Quaternion<T>(0,1,0,0)*sq);
	}
	return log(a+sqrt(a-T(1))*sqrt(a+T(1)));
}
template<typename T>
MAYBE_CONSTEXPR Quaternion<T> atanh(const Quaternion<T> &a){
	if(isnan(a)){
		return Quaternion<T>();
	}
	if(isidentity(a)){
		throw std::domain_error("atanh:singular point 1");
	}
	T w2=a.r*a.r,v2=abs(pure(a));
	if(std::abs(w2-v2-T(1))<Constants::eps<T>&&v2>Constants::eps<T>){
		throw std::domain_error("atanh:input in pathological domain");
	}
	return 0.5*log((T(1)+a)/(T(1)-a)); 
}
//string->quaternion
template<typename T>
Quaternion<T> stoq(const std::string &s){
	std::stringstream ss(s);
	Quaternion<T> q;
	if(!(ss>>q)){
		throw std::runtime_error("Failed to parse quaternion: "+s);
	}
	return q;
}
//quaternion->string
template<typename T>
std::string to_string(Quaternion<T> q){
	std::stringstream ss;
	ss<<q;
	return ss.str();
}
enum class OutMode{
	Alg,
	Tup,
	Ang
};
inline int format_storage_index(){
	static int idx=std::ios_base::xalloc();
	return idx;
}
inline std::ios_base& quat_algebraic(std::ios_base &io){
	io.iword(format_storage_index())=static_cast<long>(OutMode::Alg);
	return io;
}
inline std::ios_base& quat_tuple(std::ios_base &io){
	io.iword(format_storage_index())=static_cast<long>(OutMode::Tup);
	return io;
}
inline std::ios_base& quat_axisangle(std::ios_base &io){
	io.iword(format_storage_index())=static_cast<long>(OutMode::Ang);
	return io;
}
//a+bi+cj+dk or (a,b,c)<d
template<typename T>
std::istream &operator >>(std::istream &in,Quaternion<T> &q){
	q=Quaternion<T>();
	in>>std::ws;
	if(!in){
		return in;
	}
	std::string s;
	char c;
	while(in.get(c)){
		if(c=='\n'||c=='\r'){
			in.unget();
			break;
		}
		s.push_back(c);
	}
	s.erase(remove_if(s.begin(),s.end(),[](unsigned char c){return std::isspace(c)||c=='\n'||c=='\r';}),s.end());
	if(s.empty()){
		q=Quaternion<T>();
		return in;
	}
	std::string t="1234567890.()<ijk+-eE,";
	for(int i=0;i<s.size();i++){
		if(t.find(s[i])==t.npos){
			throw std::invalid_argument("There are illegal characters in the input");
		}
	}
	if(std::count(s.begin(),s.end(),'<')==1&&s.find("<")==s.find(")")+1&&s.find(")")!=-1&&(s.find("i")==-1&&s.find("j")==-1&&s.find("k")==-1)){
		if(s.find("(")!=0){
			s.erase(s.find("("),1);
			s.insert(0,"(");
		}
		std::string s1=s.substr(0,s.find("<"));
		std::string s2=s.substr(s.find("<")+1,s.size()-s.find("<")-1);
		s1.erase(0,1);
		s1.erase(s1.size()-1,1);
		std::string s3=s1.substr(0,s1.find(","));
		s1.erase(0,s3.size()+1);
		std::string s4=s1.substr(0,s1.find(","));
		s1.erase(0,s4.size()+1);
		Quaternion<T> u(0.0,stot<T>(s3),stot<T>(s4),stot<T>(s1));
		T theta=stot<T>(s2);
		q=std::cos(theta/2)+u*std::sin(theta/2);
		return in;
	}
	while(s.find("<")!=-1){
		s.erase(s.find("<"),1);
	}
	if(s[0]!='+'&&s[0]!='-'){
		s='+'+s;
	}
	if(s[s.size()-1]=='+'||s[s.size()-1]=='-'){
		s.erase(s.size()-1,1);
	}
	unsigned long long w=s.find("+-"),x=s.find("-+"),y=s.find("++"),z=s.find("--");
	while(w!=std::string::npos||x!=std::string::npos||y!=std::string::npos||z!=std::string::npos){
		if(w!=std::string::npos){
			s.erase(w,1);
		}
		if(x!=std::string::npos){
			s.erase(x,1);
		}
		if(y!=std::string::npos){
			s.erase(y,1);
		}
		if(z!=std::string::npos){
			s.erase(z,1);
		}
		w=s.find("+-"),x=s.find("-+"),y=s.find("++"),z=s.find("--");
	}
	if(s.empty()||(s.size()==1&&(s[0]=='+'||s[0]=='-'))){
		q=Quaternion<T>();
		return in;
	}
	w=s.find("ee"),x=s.find("eE"),y=s.find("Ee"),z=s.find("EE");
	while(w!=std::string::npos||x!=std::string::npos||y!=std::string::npos||z!=std::string::npos){
		if(w!=std::string::npos){
			s.erase(w,1);
		}
		if(x!=std::string::npos){
			s.erase(x,1);
		}
		if(y!=std::string::npos){
			s.erase(y,1);
		}
		if(z!=std::string::npos){
			s.erase(z,1);
		}
		w=s.find("ee"),x=s.find("eE"),y=s.find("Ee"),z=s.find("EE");
	}
	while(s.find("..")!=s.npos){
		s.erase(s.find(".."),1);
	}
	for(int i=0;i<s.size();i++){
		if(s[i]=='e'&&!isdigit(s[i-1])){
			s.insert(i,"0");
		}
		if(s[i]=='e'&&(!isdigit(s[i+1])&&s[i+1]!='-'&&s[i+1]!='+')){
			s.insert(i+1,"0");
		}
	}
	for(int i=0;i<s.size();i++){
		if(s[i]=='.'&&!isdigit(s[i-1])){
			s.insert(i,"0");
		}
		if(s[i]=='.'&&!isdigit(s[i+1])){
			s.insert(i+1,"0");
		}
	}
	for(int i=0;i<s.size();i++){
		if(s[i]=='E'&&!isdigit(s[i-1])){
			s.insert(i,"0");
		}
		if(s[i]=='E'&&(!isdigit(s[i+1])&&s[i+1]!='-'&&s[i+1]!='+')){
			s.insert(i+1,"0");
		}
	}
	std::vector<std::string>v;
	int start=0;
	for(int i=1;i<s.size();i++){
		if((s[i]=='+'&&s[i-1]!='e'&&s[i-1]!='E')||(s[i]=='-'&&s[i-1]!='e'&&s[i-1]!='E')){
			v.push_back(s.substr(start,i-start));
			start=i;
		}
	}
	v.push_back(s.substr(start));
	for(const auto& t:v){
		if(t.empty()){
			continue;
		}
		char sign=t[0];
		std::string coe=t.substr(1);
		if(coe.find("i")!=-1){
			std::string nn=coe.substr(0,coe.size()-1);
			T val;
			try{
				val=nn.empty()?T(1):stot<T>(nn);
			}catch(const std::exception& e){
				in.setstate(std::ios::failbit);
				in.clear();
				q=Quaternion<T>();
				throw std::invalid_argument("There are illegal characters in the input");
			}
			q.i+=(sign=='+')?val:-val;
		}else{
			if(coe.find("j")!=-1){
				std::string nn=coe.substr(0,coe.size()-1);
				T val;
				try{
					val=nn.empty()?T(1):stot<T>(nn);
				}catch(const std::exception& e){
					in.setstate(std::ios::failbit);
					in.clear();
					q=Quaternion<T>();
					throw std::invalid_argument("There are illegal characters in the input");
				}
				q.j+=(sign=='+')?val:-val;
			}else{
				if(coe.find("k")!=-1){
					std::string nn=coe.substr(0,coe.size()-1);
					T val;
					try{
						val=nn.empty()?T(1):stot<T>(nn);
					}catch(const std::exception& e){
						in.setstate(std::ios::failbit);
						in.clear();
						q=Quaternion<T>();
						throw std::invalid_argument("There are illegal characters in the input");
					}
					q.k+=(sign=='+')?val:-val;
				}else{
					T val=stot<T>(coe);
					q.r+=(sign=='+')?val:-val;
				}
			}
		}
	}
	in.clear();
	return in;
}
//a+bi+cj+dk
template<typename T>
std::ostream &operator <<(std::ostream &o,const Quaternion<T> &q){
	bool flag=0;
	auto mode=static_cast<OutMode>(o.iword(format_storage_index()));
	switch(mode){
		case OutMode::Tup:
			o<<'('<<q.r<<','<<q.i<<','<<q.j<<','<<q.k<<')';
			break;
		case OutMode::Ang:
			if(std::abs(q.r)-T(1)<=Constants::eps<T>){
				T theta=std::acos(q.r)*2;
				Quaternion<T> a=pure(q)/std::sin(theta/T(2));
				o<<'('<<a.i<<','<<a.j<<','<<a.k<<')'<<'<'<<theta;
			}
			break;
		case OutMode::Alg:default:
			if(iszero(q)){
				o<<"0";
			}else{
				if(std::abs(q.r)>Constants::eps<T>){
					o<<q.r;
					flag=1;
				}
				if(std::isnan(q.i)){
					o<<"+nani";
				}else{
					if(std::abs(q.i)>Constants::eps<T>){
						if((!flag)||q.i<T(0)){
							flag=1;
						}else{
							o<<"+";
						}
						if(std::abs(q.i-T(1))<Constants::eps<T>){
							o<<"i";
						}else{
							if(std::abs(q.i+T(1))<Constants::eps<T>){
								o<<"-i";
							}else{
								o<<q.i<<"i";
							}
						}
					}
				}
				if(std::isnan(q.j)){
					o<<"+nanj";
				}else{
					if(std::abs(q.j)>Constants::eps<T>){
						if((!flag)||q.j<T(0)){
							flag=1;
						}else{
							o<<"+";
						}
						if(std::abs(q.j-T(1))<Constants::eps<T>){
							o<<"j";
						}else{
							if(std::abs(q.j+T(1))<Constants::eps<T>){
								o<<"-j";
							}else{
								o<<q.j<<"j";
							}
						}
					}
				}
				if(std::isnan(q.k)){
					o<<"+nank";
				}else{
					if(std::abs(q.k)>Constants::eps<T>){
						if((!flag)||q.k<T(0)){
							flag=1;
						}else{
							o<<"+";
						}
						if(std::abs(q.k-T(1))<Constants::eps<T>){
							o<<"k";
						}else{
							if(std::abs(q.k+T(1))<Constants::eps<T>){
								o<<"-k";
							}else{
								o<<q.k<<"k";
							}
						}
					}
				}
			}
			break;
	}
	return o;
}
//literals
#if __cplusplus>201103L
namespace QLiterals{
	MAYBE_CONSTEXPR Quaternion<float> operator""_if(unsigned long long x){
		return Quaternion<float>((float)(0),(float)x,(float)(0),(float)(0));
	}
	MAYBE_CONSTEXPR Quaternion<float> operator""_if(long double x){
		return Quaternion<float>((float)(0),(float)x,(float)(0),(float)(0));
	}
	MAYBE_CONSTEXPR Quaternion<float> operator""_jf(unsigned long long x){
		return Quaternion<float>((float)(0),(float)(0),(float)x,(float)(0));
	}
	MAYBE_CONSTEXPR Quaternion<float> operator""_jf(long double x){
		return Quaternion<float>((float)(0),(float)(0),(float)x,(float)(0));
	}
	MAYBE_CONSTEXPR Quaternion<float> operator""_kf(unsigned long long x){
		return Quaternion<float>((float)(0),(float)(0),(float)(0),(float)x);
	}
	MAYBE_CONSTEXPR Quaternion<float> operator""_kf(long double x){
		return Quaternion<float>((float)(0),(float)(0),(float)(0),(float)x);
	}
	MAYBE_CONSTEXPR Quaternion<float> operator""_rf(unsigned long long x){
		return Quaternion<float>((float)x,(float)(0),(float)(0),(float)(0));
	}
	MAYBE_CONSTEXPR Quaternion<float> operator""_rf(long double x){
		return Quaternion<float>((float)x,(float)(0),(float)(0),(float)(0));
	}
	Quaternion<float> operator""_qf(const char* str,size_t len){
		std::string s(str,len);
		std::stringstream ss(s);
		Quaternion<float> c;
		ss>>c;
		return c;
	}
	Quaternion<float> operator""_pf(const char* str,size_t len){
		std::string s(str,len);
		if(s.find("<")==s.npos){
			return Quaternion<float>(std::stof(s));
		}
		std::string s1=s.substr(0,s.find("<"));
		std::string s2=s.substr(s.find("<")+1,s.size()-s.find("<")-1);
		s1.erase(0,1);
		s1.erase(s1.size()-1,1);
		std::string s3=s1.substr(0,s1.find(","));
		s1.erase(0,s3.size()+1);
		std::string s4=s1.substr(0,s1.find(","));
		s1.erase(0,s4.size()+1);
		Quaternion<float> u((float)(0),std::stof(s3),std::stof(s4),std::stof(s1));
		float theta=std::stof(s2);
		Quaternion<float> q=std::cos(theta/2)+u*std::sin(theta/2);
		return q;
	}
	MAYBE_CONSTEXPR Quaternion<double> operator""_id(unsigned long long x){
		return Quaternion<double>((double)(0),(double)x,(double)(0),(double)(0));
	}
	MAYBE_CONSTEXPR Quaternion<double> operator""_id(long double x){
		return Quaternion<double>((double)(0),(double)x,(double)(0),(double)(0));
	}
	MAYBE_CONSTEXPR Quaternion<double> operator""_jd(unsigned long long x){
		return Quaternion<double>((double)(0),(double)(0),(double)x,(double)(0));
	}
	MAYBE_CONSTEXPR Quaternion<double> operator""_jd(long double x){
		return Quaternion<double>((double)(0),(double)(0),(double)x,(double)(0));
	}
	MAYBE_CONSTEXPR Quaternion<double> operator""_kd(unsigned long long x){
		return Quaternion<double>((double)(0),(double)(0),(double)(0),(double)x);
	}
	MAYBE_CONSTEXPR Quaternion<double> operator""_kd(long double x){
		return Quaternion<double>((double)(0),(double)(0),(double)(0),(double)x);
	}
	MAYBE_CONSTEXPR Quaternion<double> operator""_rd(unsigned long long x){
		return Quaternion<double>((double)x,(double)(0),(double)(0),(double)(0));
	}
	MAYBE_CONSTEXPR Quaternion<double> operator""_rd(long double x){
		return Quaternion<double>((double)x,(double)(0),(double)(0),(double)(0));
	}
	Quaternion<double> operator""_qd(const char* str,size_t len){
		std::string s(str,len);
		std::stringstream ss(s);
		Quaternion<double> c;
		ss>>c;
		return c;
	}
	Quaternion<double> operator""_pd(const char* str,size_t len){
		std::string s(str,len);
		if(s.find("<")==s.npos){
			return Quaternion<double>(std::stod(s));
		}
		std::string s1=s.substr(0,s.find("<"));
		std::string s2=s.substr(s.find("<")+1,s.size()-s.find("<")-1);
		s1.erase(0,1);
		s1.erase(s1.size()-1,1);
		std::string s3=s1.substr(0,s1.find(","));
		s1.erase(0,s3.size()+1);
		std::string s4=s1.substr(0,s1.find(","));
		s1.erase(0,s4.size()+1);
		Quaternion<double> u((double)(0),std::stod(s3),std::stod(s4),std::stod(s1));
		double theta=std::stod(s2);
		Quaternion<double> q=std::cos(theta/2)+u*std::sin(theta/2);
		return q;
	}
	MAYBE_CONSTEXPR Quaternion<long double> operator""_ild(unsigned long long x){
		return Quaternion<long double>((long double)(0),(long double)x,(long double)(0),(long double)(0));
	}
	MAYBE_CONSTEXPR Quaternion<long double> operator""_ild(long double x){
		return Quaternion<long double>((long double)(0),(long double)x,(long double)(0),(long double)(0));
	}
	MAYBE_CONSTEXPR Quaternion<long double> operator""_jld(unsigned long long x){
		return Quaternion<long double>((long double)(0),(long double)(0),(long double)x,(long double)(0));
	}
	MAYBE_CONSTEXPR Quaternion<long double> operator""_jld(long double x){
		return Quaternion<long double>((long double)(0),(long double)(0),(long double)x,(long double)(0));
	}
	MAYBE_CONSTEXPR Quaternion<long double> operator""_kld(unsigned long long x){
		return Quaternion<long double>((long double)(0),(long double)(0),(long double)(0),(long double)x);
	}
	MAYBE_CONSTEXPR Quaternion<long double> operator""_kld(long double x){
		return Quaternion<long double>((long double)(0),(long double)(0),(long double)(0),(long double)x);
	}
	MAYBE_CONSTEXPR Quaternion<long double> operator""_rld(unsigned long long x){
		return Quaternion<long double>((long double)x,(long double)(0),(long double)(0),(long double)(0));
	}
	MAYBE_CONSTEXPR Quaternion<long double> operator""_rld(long double x){
		return Quaternion<long double>((long double)x,(long double)(0),(long double)(0),(long double)(0));
	}
	Quaternion<long double> operator""_qld(const char* str,size_t len){
		std::string s(str,len);
		std::stringstream ss(s);
		Quaternion<long double> c;
		ss>>c;
		return c;
	}
	Quaternion<long double> operator""_pld(const char* str,size_t len){
		std::string s(str,len);
		if(s.find("<")==s.npos){
			return Quaternion<long double>(std::stold(s));
		}
		std::string s1=s.substr(0,s.find("<"));
		std::string s2=s.substr(s.find("<")+1,s.size()-s.find("<")-1);
		s1.erase(0,1);
		s1.erase(s1.size()-1,1);
		std::string s3=s1.substr(0,s1.find(","));
		s1.erase(0,s3.size()+1);
		std::string s4=s1.substr(0,s1.find(","));
		s1.erase(0,s4.size()+1);
		Quaternion<long double> u((long double)(0),std::stold(s3),std::stold(s4),std::stold(s1));
		long double theta=std::stold(s2);
		Quaternion<long double> q=std::cos(theta/2)+u*std::sin(theta/2);
		return q;
	}
}
#endif//C++14
#endif/*QUATERNION_HPP*/
