#include <chrono>
#include <iostream>
#include "ND_array.hpp"

using namespace std;

int main()
{
    
    const int dim1 = 100;
    const int dim2 = 50;
    const int dim3 = 10;
    const int dim4 = 5;

    Array<dim1,dim2,dim3,dim4> A;

    A.fill(1.0);

    cout<<A(50,25,5,2)<<endl;
    cout<<endl<<endl<<endl;

    Array<dim1,dim2,dim3,dim4> B;
    B.fill(10.0);

    cout<<(A+B)(50,25,5,2)<<endl;
    cout<<(A+10)(50,25,5,2)<<endl;
    cout<<(A-B)(50,25,5,2)<<endl;
    cout<<(A-10)(50,25,5,2)<<endl;
    cout<<(A*B)(50,25,5,2)<<endl;
    cout<<(A*10)(50,25,5,2)<<endl;
    cout<<(A/B)(50,25,5,2)<<endl;
    cout<<(A/10)(50,25,5,2)<<endl;
    cout<<(A+=B)(50,25,5,2)<<endl;
    cout<<(A+=10)(50,25,5,2)<<endl;
    cout<<(A-=B)(50,25,5,2)<<endl;
    cout<<(A-=10)(50,25,5,2)<<endl;
    cout<<(A*=B)(50,25,5,2)<<endl;
    cout<<(A*=10)(50,25,5,2)<<endl;
    cout<<(A/=B)(50,25,5,2)<<endl;
    cout<<(A/=10)(50,25,5,2)<<endl;
    
    cout<<endl<<endl<<endl;

    Array<dim1,dim2,dim3,dim4> C;
    C = -B;
    cout<<C(50,25,5,2)<<endl;
    cout<<-B(50,25,5,2)<<endl;
    cout<<endl<<endl<<endl;

    cout<<(A-B).abs()(50,25,5,2)<<endl;
    cout<<abs(A-B)(50,25,5,2)<<endl;
    cout<<endl<<endl<<endl;
    
    cout<<(A-B).min()<<endl;
    cout<<min(A-B)<<endl;
    cout<<(A-B).max()<<endl;
    cout<<max(A-B)<<endl;
    cout<<endl<<endl<<endl;
    

    /*// 1D array
    const int dim = 100;

    Array<dim> A;

    A.fill(1.0);

    cout<<A(50)<<endl;
    cout<<endl<<endl<<endl;

    Array<dim> B;
    B.fill(10.0);

    cout<<(A+B)(50)<<endl;
    cout<<(A+10)(50)<<endl;
    cout<<(A-B)(50)<<endl;
    cout<<(A-10)(50)<<endl;
    cout<<(A*B)(50)<<endl;
    cout<<(A*10)(50)<<endl;
    cout<<(A/B)(50)<<endl;
    cout<<(A/10)(50)<<endl;
    cout<<(A+=B)(50)<<endl;
    cout<<(A+=10)(50)<<endl;
    cout<<(A-=B)(50)<<endl;
    cout<<(A-=10)(50)<<endl;
    cout<<(A*=B)(50)<<endl;
    cout<<(A*=10)(50)<<endl;
    cout<<(A/=B)(50)<<endl;
    cout<<(A/=10)(50)<<endl;
    
    cout<<endl<<endl<<endl;

    Array<dim> C;
    C = -B;
    cout<<C(50)<<endl;
    cout<<-B(50)<<endl;
    cout<<endl<<endl<<endl;

    cout<<(A-B).abs()(50)<<endl;
    cout<<abs(A-B)(50)<<endl;
    cout<<endl<<endl<<endl;
    
    cout<<(A-B).min()<<endl;
    cout<<min(A-B)<<endl;
    cout<<(A-B).max()<<endl;
    cout<<max(A-B)<<endl;
    cout<<endl<<endl<<endl;
    */
    return 0;
}