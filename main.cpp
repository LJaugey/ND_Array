#include <chrono>
#include <iostream>
#include "ND_array.hpp"

using namespace std;

int main()
{
    Array<5,4,3,2> A;
    Array<5,4,3,2> B;
    
    cout<<A[1][2]<<endl;
    A.fill(1.0);
    cout<<A[1][2]<<endl;
    
    B.fill(3.0);


    Array<5,4,3,2> C = A+B;

    C(1,2,0,0) = 42.0;

    cout<<C[1][2]<<endl;

    return 0;
}