#include "Quaternion.hpp"
#include <cassert>
Quaternion<double> a,b,c,d,e;
Quaternion<long double> x,y,u,v,w,cc,cd,ce;
int main(){
  cin>>a>>b;
  x=a;y=b;
  c=a+b;d=a-b;e=a*b;
  u=x+y;v=x-y;w=x*y;
  cc=c;cd=d;ce=e;
  assert(almost_equal(cc,u,1e-6));
  assert(almost_equal(cd,v,1e-6));
  assert(almost_equal(ce,w,1e-6));
}
