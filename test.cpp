#include <bits/stdc++.h>
#include <cassert>
#include "Quaternion.hpp"
using namespace std;
int main(){
    //1. String parsing test
    auto q=stoq<double>(1+2i+3j+4k);
    assert(q.real()==1.0);
    assert(q.imagi()==2.0);
    assert(q.imagj()==3.0);
    assert(q.imagk()==4.0);
    cout<<"String parsing test passed\n";
    //2. Rotation test
    auto rot=Quaternion<double>::FromAxisAngle(1.5708,0,0,1);
    auto v=rot.RotateVector(1,0,0);
    assert(abs(v[0])<1e-6&&abs(v[1]-1)<1e-6);
    cout<<"Rotation test passed\n";
    //3.Output format test
    cout<<quat_tuple<<q<<endl;
    cout<<quat_axisangle<<q<<endl;
    cout<<q<<endl;
    cout<<"All tests passed!";
    return 0;
}
