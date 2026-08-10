#include "Quaternion.hpp"
#include <iostream>
using namespace QLiterals;
int main(){
  auto direc=std::array<float,3>{1.0f,0.0f,0.0f};
  auto q=Quaternion<float>::LookAt(direc);
  std::cout<<"LookAt(1,0,1)="<<q<<std::endl;
  return 0;
}
