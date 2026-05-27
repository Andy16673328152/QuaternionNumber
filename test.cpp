#include <iostream>
#include <cassert>
#include "Quaternion.hpp"
int main(){
    //1. String parsing test
    auto q=stoq<double>("1+2i+3j+4k");
    assert(q.real()==1.0);
    assert(q.imagi()==2.0);
    assert(q.imagj()==3.0);
    assert(q.imagk()==4.0);
    std::cout<<"String parsing test passed\n";
    //2. Rotation test
    auto rot=Quaternion<double>::FromAxisAngle(1.5708,0,0,1);
    auto v=rot.RotateVector(1,0,0);
    assert(std::abs(v[0])<1e-6&&std::abs(v[1]-1)<1e-6);
    std::cout<<"Rotation test passed\n";
    //3.Output format test
    std::cout<<quat_tuple<<q<<"\n";
    std::cout<<quat_axisangle<<q<<"\n";
    std::cout<<q<<"\n";
    std::cout<<"All tests passed!";
    return 0;
}
