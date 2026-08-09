# QuaternionNumber
This library is a header-only library,it provides:

    - Quaternion basic operations;
    
    - Transcendental functions;
    
    - Literals;
    
    - `slerp` and `nlerp` interpolation(SIMD);

    - `squad` interpolation;

    - distances;

    - Conversion functions with other libraries.
    

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![C++](https://img.shields.io/badge/C%2B%2B-17/20/23-blue.svg)](https://isocpp.org/)
[![Version](https://img.shields.io/badge/2.1.0-green.svg)]()
[![Header-Only](https://img.shields.io/badge/header-only-green.svg)]()

Example code:
```cpp
//You can copy this code.
#include "Quaternion.hpp"
#include <iostream>
using namespace QLiterals;
int main(){
    auto q=stoq<double>("1+2i+3j+4k");
    std::cout<<q<<std::endl;//1+2i+3j+4k
}
```
