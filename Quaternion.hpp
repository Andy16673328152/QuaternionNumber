#ifndef QUATERNION_HPP
#define QUATERNION_HPP
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
namespace Constants{
	template<typename T>
	inline constexpr T PI=T(3.14159265358979323846);
	template<typename T>
	inline constexpr T eps=T(1e-10);
}
std::string TrimTrailingZeros(std::string a){
	if(a.empty()){
		return a;
	}
        while(a[a.size()-1]=='0'){
            a.erase(a.size()-1,1);
        }
    return a;
}
template <typename T>
struct Quaternion{
	static_assert(std::is_same<T,float>::value||std::is_same<T,double>::value||std::is_same<T,long double>::value,"T must be 'float','double' or 'long double',not 'int'");
	T r,i,j,k;
	constexpr Quaternion conj()const{
		return {r,-i,-j,-k};
	}
	constexpr T real()const{
		return r;
	}
	void real(T val){
		r=val;
	}
	constexpr T imagi()const{
		return i;
	}
	void imagi(T val){
		i=val;
	}
	constexpr T imagj()const{
		return j;
	}
	void imagj(T val){
		j=val;
	}
	constexpr T imagk()const{
		return k;
	}
	void imagk(T val){
		k=val;
	}
	constexpr Quaternion():r(0.0),i(0.0),j(0.0),k(0.0){}
	constexpr Quaternion(T a,T b,T c,T d):r(a),i(b),j(c),k(d){}
	constexpr Quaternion(T a):r(a){}
	explicit constexpr Quaternion(const std::complex<T> &a):r(a.real()),i(a.imag()),j(0),k(0){}
	constexpr Quaternion(const std::array<T,4> &a):r(a[0]),i(a[1]),j(a[2]),k(a[3]){}
	constexpr Quaternion(const std::array<T,3> &a):r(0),i(a[0]),j(a[1]),k(a[2]){}
	constexpr Quaternion(const std::initializer_list<T> &l){
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
	Quaternion(const Quaternion<U> &a):r(a.r),i(a.i),j(a.j),k(a.k){}
	Quaternion operator =(Quaternion a){
		r=a.r;
		i=a.i;
		j=a.j;
		k=a.k;
		return *this;
	}
	template <typename U>
	Quaternion operator =(const Quaternion<U> &a){
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
		return Quaternion(1.0);
	}
	static Quaternion unit(T r,T i,T j,T k){
		Quaternion q(r,i,j,k);
		return q.normalized();
	}
	Quaternion& normalize(){
		T al=std::sqrt(r*r+i*i+j*j+k*k);
		if(al>Constants::eps<T>){
			r/=al;
			i/=al;
			j/=al;
			k/=al;
		}
		return *this;
	}
	Quaternion normalized()const{
		Quaternion re(*this);
		return re.normalize();
	}
	bool isunit(T eeps=Constants::eps<T>)const{
		T alsq=r*r+i*i+j*j+k*k;
		return std::abs(alsq-1)<eeps;
	}
	static Quaternion FromAxisAngle(T angle,T x,T y,T z){
		Quaternion q(std::cos(angle/2),x*std::sin(angle/2),y*std::sin(angle/2),z*std::sin(angle/2));
		return q.normalized();
	} 
	std::array<T,9> ToRotationMatrix()const{
		T xx=r*r,yy=i*i,zz=j*j,ww=k*k,xy=r*i,xz=r*j,xw=r*k,yz=i*j,yw=i*k,zw=j*k;
		return {
			1-2*(yy+zz),2*(xy-zw),2*(xz+yw),
			2*(xy+zw),1-2*(xx+zz),2*(yz-xw),
			2*(xz-yw),2*(yz+xw),1-2*(xx+yy)
		};
	}
	static Quaternion FromRotationMatrix(const std::array<T,9>& mat){
		T t=mat[0]+mat[4]+mat[8];
		Quaternion q;
		if(t>0){
	 		T s=0.5*(std::sqrt(t+1));
	 		q.r=0.25/s;
	 		q.i=(mat[5]-mat[7])*s;
	 		q.j=(mat[6]-mat[2])*s;
	 		q.k=(mat[1]-mat[3])*s;
		}else{
		 	if(mat[0]>mat[4]&&mat[0]>mat[8]){
		 		T s=2*std::sqrt(1+mat[0]-mat[4]-mat[8]);
		 		q.r=(mat[5]-mat[7])/s;
		 		q.i=0.25*s;
		 		q.j=(mat[3]+mat[1])/s;
		 		q.k=(mat[6]+mat[2])/s;
			}else{
				if(mat[4]>mat[8]){
					T s=2*std::sqrt(1+mat[4]-mat[0]-mat[8]);
					q.r=(mat[6]-mat[2])/s;
					q.i=(mat[3]+mat[1])/s;
					q.j=0.25*s;
					q.k=(mat[7]+mat[5])/s;
				}else{
					T s=2*std::sqrt(1+mat[8]-mat[4]-mat[0]);
					q.r=(mat[1]-mat[3])/s;
					q.i=(mat[6]+mat[2])/s;
					q.j=(mat[7]+mat[5])/s;
					q.k=0.25*s;
				}
			}
		}
		return q;
	}
	enum class EulerOrder{
		XYZ,XZY,YXZ,YZX,ZXY,ZYX
	};
	std::array<T,3> ToEulerAngles(EulerOrder order=EulerOrder::XYZ)const{
		Quaternion q=this->normalized();
		T w=q.r,x=q.i,y=q.j,z=q.k;
		T roll,pitch,yaw;
		const T epss=T(0.999999);
		auto clasin=[](T val){
			return std::asin(std::max(T(-1),std::min(T(1),val)));
		};
		switch(order){
			case EulerOrder::XYZ:
				pitch=clasin(2*(w*y-z*x));
				if(std::abs(std::sin(pitch))<epss){
					roll=atan2(2*(w*x+y*z),1-2*(x*x+y*y));
					yaw=atan2(2*(w*z+x*y),1-2*(y*y+z*z));
				}else{
					roll=atan2(2*(x*z-w*y),1-2*(x*x+z*z));
					yaw=0;
				}
				return {roll,pitch,yaw};
			case EulerOrder::XZY:
				yaw=clasin(2*(w*z+x*y));
				if(std::abs(std::sin(yaw))<epss){
					pitch=atan2(2*(w*y-x*z),1-2*(y*y+z*z));
					roll=atan2(2*(w*x-y*z),1-2*(x*x+z*z));
				}else{
					roll=atan2(2*(x*z+w*y),1-2*(x*x+y*y));
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
					roll=atan2(2*(w*x+y*z),1-2*(x*x+y*y));
					yaw=atan2(2*(w*z+x*y),1-2*(y*y+z*z));
				}else{
					yaw=atan2(2*(x*z-w*y),1-2*(x*x+z*z));
					roll=0;
				}
				return {roll,pitch,yaw};
			default:
				return {0,0,0};
		}
	}
	static Quaternion FromEulerAngles(T roll,T pitch,T yaw,EulerOrder order=EulerOrder::XYZ){
		T c1=std::cos(roll/2);
		T s1=std::sin(roll/2);
		T c2=std::cos(pitch/2);
		T s2=std::sin(pitch/2);
		T c3=std::cos(yaw/2);
		T s3=std::sin(yaw/2);
		switch(order){
			case EulerOrder::XYZ:
				Quaternion q(
					c1*c2*c3-s1*s2*s3,
					s1*c2*c3+c1*s2*s3,
					c1*s2*c3-s1*c2*s3,
					c1*c2*s3+s1*s2*c3
				);
				return q.normalized();
			case EulerOrder::XZY:
				Quaternion q1(
					c1*c2*c3+s1*s2*s3,
					s1*c2*c3-c1*s2*s3,
					c1*c2*s3-s1*s2*c3,
					c1*s2*c3+s1*c2*s3
				);
				return q1.normalized();
			case EulerOrder::YXZ:
				Quaternion q2(
					c1*c2*c3+s1*s2*s3,
					s1*c2*c3+c1*s2*s3,
					c1*s2*c3-s1*c2*s3,
					c1*c2*s3-s1*s2*c3
				);
				return q2.normalized();
			case EulerOrder::YZX:
				Quaternion q3(
					c1*c2*c3-s1*s2*s3,
					c1*s2*c3+s1*c2*s3,
					s1*c2*c3+c1*s2*s3,
					c1*c2*s3-s1*s2*c3
				);
				return q3.normalized();
			case EulerOrder::ZXY:
				Quaternion q4(
					c1*c2*c3-s1*s2*s3,
					c1*s2*c3-s1*c2*s3,
					c1*c2*s3+s1*s2*c3,
					s1*c2*c3+c1*s2*s3
				);
				return q4.normalized();
			case EulerOrder::ZYX:
				Quaternion q5(
					c1*c2*c3+s1*s2*s3,
					c1*c2*s3-s1*s2*c3,
					c1*s2*c3+s1*c2*s3,
					s1*c2*c3-c1*s2*s3
				);
				return q5.normalized();
			default:
				return Quaternion(1,0,0,0);
		}
	}
	std::array<T,3> RotateVector(T x,T y,T z)const{
		T rx=2*(j*z-k*y),ry=2*(k*x-i*z),rz=2*(i*y-j*x);
		return {x+r*rx+(j*rz-k*ry),y+r*ry+(k*rx-i*rz),z+r*rz+(i*ry-j*rx)};
	}
	std::array<T,3> RotateVector(const std::array<T,3> &v)const{
		return RotateVector(v[0],v[1],v[2]);
	}
	std::array<T,3> axis()const{
		if(std::isnan(r)||std::isnan(i)||std::isnan(j)||std::isnan(k)){
			return {1,0,0};
		}
		T cr=std::max(T(-1),std::min(T(1),r));
		T a=std::sqrt(1-cr*cr);
		if(a<Constants::eps<T>||std::isnan(a)){
			return {1,0,0};
		}
		return {i/a,j/a,k/a};
	}
	T angle()const{
		if(std::isnan(r)||std::isnan(i)||std::isnan(j)||std::isnan(k)){
			return T(0);
		}
		T cr=std::max(T(-1),std::min(T(1),r));
		return 2*std::acos(cr);
	}
	T angleDegrees()const{
		return angle()*180/Constants::PI<T>;
	}
	Quaternion inverse()const{
		return (this->conj())/(r*r+i*i+j*j+k*k);
	}
};
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
template<typename T>
Quaternion<T> operator +(Quaternion<T> a,T b){
	a.r+=b;
	return a;
} 
template<typename T>
Quaternion<T> operator +(T a,Quaternion<T> b){
	b.r+=a;
	return b;
}
template<typename T>
Quaternion<T> operator -(Quaternion<T> a,T b){
	a.r-=b;
	return a;
} 
template<typename T>
Quaternion<T> operator -(Quaternion<T> a){
	return {-a.r,-a.i,-a.j,-a.k};
}
template<typename T>
Quaternion<T> operator -(T a,Quaternion<T> b){
	b.r-=a;
	b=-b;
	return b;
}
template<typename T>
Quaternion<T> operator *(Quaternion<T> a,T b){
	a.r*=b;
	a.i*=b;
	a.j*=b;
	a.k*=b;
	return a;
}
template<typename T>
Quaternion<T> operator *(T a,Quaternion<T> b){
	b.r*=a;
	b.i*=a;
	b.j*=a;
	b.k*=a;
	return b;
}
template<typename T>
Quaternion<T> operator /(Quaternion<T> a,T b){
	if(std::abs(b)<Constants::eps<T>){
		throw std::runtime_error("Division by zero");
	}
	a.r/=b;
	a.i/=b;
	a.j/=b;
	a.k/=b;
	return a;
}
template<typename T>
T abs(Quaternion<T> a){
	return std::sqrt(a.r*a.r+a.i*a.i+a.j*a.j+a.k*a.k);
}
template<typename T>
Quaternion<T> operator /(T a,Quaternion<T> b){
	if(abs(b)<Constants::eps<T>){
		throw std::runtime_error("Division by zero");
	}
	return a*(b.conj()/(abs(b)*abs(b)));
}
template<typename T>
Quaternion<T> operator +(Quaternion<T> a,Quaternion<T> b){
	a.r+=b.r;
	a.i+=b.i;
	a.j+=b.j;
	a.k+=b.k;
	return a;
}
template<typename T>
Quaternion<T> operator -(Quaternion<T> a,Quaternion<T> b){
	a.r-=b.r;
	a.i-=b.i;
	a.j-=b.j;
	a.k-=b.k;
	return a;
}
template<typename T>
Quaternion<T> operator *(Quaternion<T> a,Quaternion<T> b){
	return Quaternion<T>(a.r*b.r-a.i*b.i-a.j*b.j-a.k*b.k,a.r*b.i+a.i*b.r+a.j*b.k-a.k*b.j,a.r*b.j+a.j*b.r+a.k*b.i-a.i*b.k,a.r*b.k+a.k*b.r+a.i*b.j-a.j*b.i);
}
template<typename T>
Quaternion<T> operator /(Quaternion<T> a,Quaternion<T> b){
	return a*(1.0/b);
}
template<typename T>
bool operator ==(Quaternion<T> a,Quaternion<T> b){
	if(std::isnan(a.r)||std::isnan(a.i)||std::isnan(a.j)||std::isnan(a.k)||std::isnan(b.r)||std::isnan(b.i)||std::isnan(b.j)||std::isnan(b.k)){
		return false;
	}
	if(std::isinf(a.r)||std::isinf(a.i)||std::isinf(a.j)||std::isinf(a.k)||std::isinf(b.r)||std::isinf(b.i)||std::isinf(b.j)||std::isinf(b.k)){
		return (a.r==b.r)&&(a.i==b.i)&&(a.j==b.j)&&(a.k==b.k);
	}
	return (std::abs(a.r-b.r)<Constants::eps<T>)&&(std::abs(a.i-b.i)<Constants::eps<T>)&&(std::abs(a.j-b.j)<Constants::eps<T>)&&(std::abs(a.k-b.k)<Constants::eps<T>);
}
template<typename T>
bool operator !=(Quaternion<T> a,Quaternion<T> b){
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
template<typename T>
Quaternion<T> pure(Quaternion<T> a){
	return Quaternion<T>(0.0,a.i,a.j,a.k);
}
template<typename T>
Quaternion<T> log(Quaternion<T> a){
	T pa=abs(pure(a));
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::log(a.r),0,0,0);
	}
	Quaternion<T> u(0.0,a.i/pa,a.j/pa,a.k/pa);
	return std::log(abs(a))+u*std::acos(a.r/abs(a));
}
template<typename T>
Quaternion<T> log10(Quaternion<T> a){
	return log(a)/std::log(10);
}
template<typename T>
Quaternion<T> log2(Quaternion<T> a){
	return log(a)/std::log(2);
}
template<typename T>
Quaternion<T> logp(Quaternion<T> a,Quaternion<T> b){
	return log(a)/log(b);
}
template<typename T>
Quaternion<T> exp(Quaternion<T> a){
	T pa=abs(pure(a));
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::exp(a.r),0,0,0);
	}
	Quaternion<T> u(0.0,a.i/pa,a.j/pa,a.k/pa);
	return std::exp(a.r)*(std::cos(pa)+u*std::sin(pa));
}
template<typename T>
Quaternion<T> pow(Quaternion<T> a,Quaternion<T> b){
	return exp(log(a)*b);
}
template<typename T>
Quaternion<T> nthrt(Quaternion<T> a,Quaternion<T> b){
	return pow(a,1.0/b);
}
template<typename T>
Quaternion<T> sqrt(Quaternion<T> a){
	return nthrt(a,Quaternion<T>(2.0));
}
template<typename T>
Quaternion<T> cbrt(Quaternion<T> a){
	return nthrt(a,Quaternion<T>(3.0));
}
template<typename T>
T dot(Quaternion<T> a,Quaternion<T> b){
	return a.r*b.r+a.i*b.i+a.j*b.j+a.k*b.k;
}

template <typename T>
Quaternion<T> unit(T b,T c,T d){
	return Quaternion<T>(0.0,b/std::sqrt(b*b+c*c+d*d),c/std::sqrt(b*b+c*c+d*d),d/std::sqrt(b*b+c*c+d*d));
}
template<typename T>
Quaternion<T> slerp(const Quaternion<T> &a,const Quaternion<T> &b,T t){
	T d=dot(a,b),sign=1;
	if(d<0){
		d=-d;
		sign=-1;
	}
	const T dt=T(0.9995);
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
template<typename T>
Quaternion<T> omega(Quaternion<T> a,int n,int k){
	k%=n;
	if(a.is_real()){
		if(std::abs(a.r)<Constants::eps<T>){
			return Quaternion(0.0);
		}else{
			if(a.r>0){
				T poww=std::pow(a.r,1.0/n);
				Quaternion<T> q(std::cos(2*Constants::PI<T>*k/n),std::sin(2*Constants::PI<T>*k/n),0,0);
				return poww*q;
			}else{
				T poww=std::pow(-a.r,1.0/n);
				Quaternion<T> q(std::cos((Constants::PI<T>+2*Constants::PI<T>*k)/n),std::sin((Constants::PI<T>+2*Constants::PI<T>*k)/n),0,0);
				return poww*q;
			}
		}
	}
	T r=abs(a),theta=a.angle();
	Quaternion<T> u=pure(a)/std::sqrt(a.i*a.i+a.j*a.j+a.k*a.k);
	return std::pow(r,1.0/n)*(std::cos((theta+2*Constants::PI<T>*k)/n)+u*std::sin((theta+2*Constants::PI<T>*k)/n));
}
template<typename T>
Quaternion<T> sin(Quaternion<T> a){
	T pa=abs(pure(a));
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::sin(a.r),0,0,0);
	}
	Quaternion<T> u(0.0,a.i/pa,a.j/pa,a.k/pa);
	return std::sin(a.r)*std::cosh(pa)+u*std::cos(a.r)*std::sinh(pa);
}
template<typename T>
Quaternion<T> cos(Quaternion<T> a){
	T pa=abs(pure(a));
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::cos(a.r),0,0,0);
	}
	Quaternion<T> u(0.0,a.i/pa,a.j/pa,a.k/pa);
	return std::cos(a.r)*std::cosh(pa)-u*std::sin(a.r)*std::sinh(pa);
}
template<typename T>
Quaternion<T> tan(Quaternion<T> a){
	return sin(a)/cos(a);
}
template<typename T>
Quaternion<T> asin(Quaternion<T> a){
	T pa=abs(pure(a));
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::asin(a.r),0,0,0);
	}
	Quaternion<T> u(0.0,a.i/pa,a.j/pa,a.k/pa);
	return -u*log(a*u+sqrt(1.0-a*a));
}
template<typename T>
Quaternion<T> acos(Quaternion<T> a){
	T pa=abs(pure(a));
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::acos(a.r),0,0,0);
	}
	Quaternion<T> u(0.0,a.i/pa,a.j/pa,a.k/pa);
	return -u*log(a+u*sqrt(1.0-a*a));
}
template<typename T>
Quaternion<T> atan(Quaternion<T> a){
	T pa=abs(pure(a));
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::atan(a.r),0,0,0);
	}
	Quaternion<T> u(0.0,a.i/pa,a.j/pa,a.k/pa);
	return u/2*log((u+a)/(u-a));
}
template<typename T>
Quaternion<T> sinh(Quaternion<T> a){
	T pa=abs(pure(a));
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::sinh(a.r),0,0,0);
	}
	Quaternion<T> u(0.0,a.i/pa,a.j/pa,a.k/pa);
	return std::sinh(a.r)*std::cos(pa)+u*std::cosh(a.r)*std::sin(pa);
}
template<typename T>
Quaternion<T> cosh(Quaternion <T> a){
	T pa=abs(pure(a));
	if(pa<Constants::eps<T>){
		return Quaternion<T>(std::cosh(a.r),0,0,0);
	}
	Quaternion<T> u(0.0,a.i/pa,a.j/pa,a.k/pa);
	return std::cosh(a.r)*std::cos(pa)+u*std::sinh(a.r)*std::sin(pa);
}
template<typename T>
Quaternion<T> tanh(Quaternion<T> a){
	return sinh(a)/cosh(a);
}
template<typename T>
Quaternion<T> asinh(Quaternion<T> a){
	return log(a+sqrt(a*a+1));
}
template<typename T>
Quaternion<T> acosh(Quaternion<T> a){
	return log(a+sqrt(a-1.0)*sqrt(a+1.0));
}
template<typename T>
Quaternion<T> atanh(Quaternion<T> a){
	return 0.5*log((1.0+a)/(1.0-a)); 
}
template<typename T>
Quaternion<T> stoq(const std::string &s){
	std::stringstream ss(s);
	Quaternion<T> q;
	if(!(ss>>q)){
		throw std::runtime_error("Failed to parse quaternion: "+s);
	}
	return q;
}
template<typename T>
std::string to_string(Quaternion<T> q){
	std::stringstream ss;
	ss<<q;
	return ss.str();
}
template<typename T>
std::istream &operator >>(std::istream &in,Quaternion<T> &q){
	std::string s;
	in>>s;
	s.erase(remove_if(s.begin(),s.end(),::isspace),s.end());
	if(s.empty()){
		q=Quaternion<T>(0.0,0.0,0.0,0.0);
		return in;
	}
	if(s.find("<")!=-1){
		std::string s1=s.substr(0,s.find("<"));
		std::string s2=s.substr(s.find("<")+1,s.size()-s.find("<")-1);
		s1.erase(0,1);
		s1.erase(s1.size()-1,1);
		std::string s3=s1.substr(0,s1.find(","));
		s1.erase(0,s3.size()+1);
		std::string s4=s1.substr(0,s1.find(","));
		s1.erase(0,s4.size()+1);
		Quaternion<T> u(0.0,std::stold(s3),std::stold(s4),std::stold(s1));
		T theta=stold(s2);
		q=std::cos(theta/2)+u*std::sin(theta/2);
		return in;
	}
	if(s[0]!='+'&&s[0]!='-'){
		s='+'+s;
	}
	if(s[s.size()-1]=='+'||s[s.size()-1]=='-'){
		s.erase(s.size()-1,1);
	}
	while(s.find("++")!=std::string::npos){
		s.erase(s.find("++"),1);
	}
	while(s.find("--")!=std::string::npos){
		s.erase(s.find("--"),1);
	}
	std::vector<std::string>v;
	int start=0;
	for(int i=1;i<s.size();i++){
		if(s[i]=='+'||(s[i]=='-'&&s[i-1]!='e'&&s[i-1]!='E')){
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
				val=nn.empty()?1.0L:std::stold(nn);
			}catch(const std::exception& e){
				in.setstate(std::ios::failbit);
				return in;
			}
			q.i+=(sign=='+')?val:-val;
		}else{
			if(coe.find("j")!=-1){
				std::string nn=coe.substr(0,coe.size()-1);
				T val;
				try{
					val=nn.empty()?1.0L:std::stold(nn);
				}catch(const std::exception& e){
					in.setstate(std::ios::failbit);
					return in;
				}
				q.j+=(sign=='+')?val:-val;
			}else{
				if(coe.find("k")!=-1){
					std::string nn=coe.substr(0,coe.size()-1);
					T val;
					try{
						val=nn.empty()?1.0L:std::stold(nn);
					}catch(const std::exception& e){
						in.setstate(std::ios::failbit);
						return in;
					}
					q.k+=(sign=='+')?val:-val;
				}else{
					T val=std::stold(coe);
					q.r+=(sign=='+')?val:-val;
				}
			}
		}
	}
	return in;
}
template<typename T>
std::ostream &operator <<(std::ostream &o,const Quaternion<T> &q){
	std::string s1=std::to_string(q.r),s2=std::to_string(q.i),s3=std::to_string(q.j),s4=std::to_string(q.k);
	s1=TrimTrailingZeros(s1);
	s2=TrimTrailingZeros(s2);
	s3=TrimTrailingZeros(s3);
	s4=TrimTrailingZeros(s4);
	bool flag=0;
	if(s1.find(".")==s1.size()-1){
		s1.erase(s1.size()-1,1);
	}
	if(s2.find(".")==s2.size()-1){
		s2.erase(s2.size()-1,1);
	}
	if(s3.find(".")==s3.size()-1){
		s3.erase(s3.size()-1,1);
	}
	if(s4.find(".")==s4.size()-1){
		s4.erase(s4.size()-1,1);
	}
	if(std::abs(q.r)<Constants::eps<T>&&std::abs(q.i)<Constants::eps<T>&&std::abs(q.j)<Constants::eps<T>&&std::abs(q.k)<Constants::eps<T>){
		o<<"0";
		return o;
	}
	if(std::abs(q.r)>Constants::eps<T>){
		o<<s1;
		flag=1;
	}
	if(std::abs(q.i)>Constants::eps<T>){
		if((!flag)||q.i<0){
			flag=1;
		}else{
			o<<"+";
		}
		if(std::abs(q.i-1)<Constants::eps<T>){
			o<<"i";
		}else{
			if(std::abs(q.i+1)<Constants::eps<T>){
				o<<"-i";
			}else{
				o<<s2<<"i";
			}
		}
	}
	if(std::abs(q.j)>Constants::eps<T>){
		if((!flag)||q.j<0){
			flag=1;
		}else{
			o<<"+";
		}
		if(std::abs(q.j-1)<Constants::eps<T>){
			o<<"j";
		}else{
			if(std::abs(q.j+1)<Constants::eps<T>){
				o<<"-j";
			}else{
				o<<s3<<"j";
			}
		}
	}
	if(std::abs(q.k)>Constants::eps<T>){
		if((!flag)||q.k<0){
			flag=1;
		}else{
			o<<"+";
		}
		if(std::abs(q.k-1)<Constants::eps<T>){
			o<<"k";
		}else{
			if(std::abs(q.k+1)<Constants::eps<T>){
				o<<"-k";
			}else{
				o<<s4<<"k";
			}
		}
	}
	return o;
}
#if __cplusplus>201103L
namespace QLiterals{
	Quaternion<float> operator""_if(unsigned long long x){
		return Quaternion<float>(0,(float)x,0,0);
	}
	Quaternion<float> operator""_if(long double x){
		return Quaternion<float>(0,(float)x,0,0);
	}
	Quaternion<float> operator""_jf(unsigned long long x){
		return Quaternion<float>(0,0,(float)x,0);
	}
	Quaternion<float> operator""_jf(long double x){
		return Quaternion<float>(0,0,(float)x,0);
	}
	Quaternion<float> operator""_kf(unsigned long long x){
		return Quaternion<float>(0,0,0,(float)x);
	}
	Quaternion<float> operator""_kf(long double x){
		return Quaternion<float>(0,0,0,(float)x);
	}
	Quaternion<float> operator""_rf(unsigned long long x){
		return Quaternion<float>((float)x,0,0,0);
	}
	Quaternion<float> operator""_rf(long double x){
		return Quaternion<float>((float)x,0,0,0);
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
		Quaternion<float> u(0.0,std::stof(s3),std::stof(s4),std::stof(s1));
		float theta=std::stof(s2);
		Quaternion<float> q=std::cos(theta/2)+u*std::sin(theta/2);
		return q;
	}
	Quaternion<double> operator""_id(unsigned long long x){
		return Quaternion<double>(0,(double)x,0,0);
	}
	Quaternion<double> operator""_id(long double x){
		return Quaternion<double>(0,(double)x,0,0);
	}
	Quaternion<double> operator""_jd(unsigned long long x){
		return Quaternion<double>(0,0,(double)x,0);
	}
	Quaternion<double> operator""_jd(long double x){
		return Quaternion<double>(0,0,(double)x,0);
	}
	Quaternion<double> operator""_kd(unsigned long long x){
		return Quaternion<double>(0,0,0,(double)x);
	}
	Quaternion<double> operator""_kd(long double x){
		return Quaternion<double>(0,0,0,(double)x);
	}
	Quaternion<double> operator""_rd(unsigned long long x){
		return Quaternion<double>((double)x,0,0,0);
	}
	Quaternion<double> operator""_rd(long double x){
		return Quaternion<double>((double)x,0,0,0);
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
		Quaternion<double> u(0.0,std::stod(s3),std::stod(s4),std::stod(s1));
		double theta=std::stod(s2);
		Quaternion<double> q=std::cos(theta/2)+u*std::sin(theta/2);
		return q;
	}
	Quaternion<long double> operator""_ild(unsigned long long x){
		return Quaternion<long double>(0,(long double)x,0,0);
	}
	Quaternion<long double> operator""_ild(long double x){
		return Quaternion<long double>(0,(long double)x,0,0);
	}
	Quaternion<long double> operator""_jld(unsigned long long x){
		return Quaternion<long double>(0,0,(long double)x,0);
	}
	Quaternion<long double> operator""_jld(long double x){
		return Quaternion<long double>(0,0,(long double)x,0);
	}
	Quaternion<long double> operator""_kld(unsigned long long x){
		return Quaternion<long double>(0,0,0,(long double)x);
	}
	Quaternion<long double> operator""_kld(long double x){
		return Quaternion<long double>(0,0,0,(long double)x);
	}
	Quaternion<long double> operator""_rld(unsigned long long x){
		return Quaternion<long double>((long double)x,0,0,0);
	}
	Quaternion<long double> operator""_rld(long double x){
		return Quaternion<long double>((long double)x,0,0,0);
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
		Quaternion<long double> u(0.0,std::stold(s3),std::stold(s4),std::stold(s1));
		long double theta=std::stold(s2);
		Quaternion<long double> q=std::cos(theta/2)+u*std::sin(theta/2);
		return q;
	}
}
#endif//C++14
#endif/*QUATERNION_HPP*/
