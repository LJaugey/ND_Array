#include <iostream>
#include "linear_algebra.hpp"


using namespace std;

int main()
{
    Matrix<4,4> A;
    A.fill(0.5);
    Matrix<4,4> B;
    B.fill(0.5);

    Matrix<4,4> C = A*B;

    cout<<A<<endl<<endl;
    cout<<B<<endl<<endl;
    cout<<C<<endl<<endl;

    cout<<endl<<endl;

    Vector<4> V;
    V.fill(0.5);

    Vector<4> res = A*V;

    cout<<res<<endl<<endl;
    cout<<endl<<endl;

    cout<<V*res<<endl;

    return 0;
}