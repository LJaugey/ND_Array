#include <chrono>
#include <iostream>
#include <type_traits>
#include "ND_array.hpp"

using namespace std;

int main()
{
    Array<2,3> X(2);
    Array<2,3> Y;

    Y = 1+X+1;

    cout<<Y<<endl;

    
    const int dim1 = 100;
    const int dim2 = 50;
    const int dim3 = 10;
    const int dim4 = 5;

    Array<dim1,dim2,dim3,dim4> A(1.0);
    Array<dim1,dim2,dim3,dim4> B(10.0);

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

    cout<<((A-B*C)/A).abs()(50,25,5,2)<<endl;
    cout<<abs((A-B*C)/A)(50,25,5,2)<<endl;
    cout<<endl<<endl<<endl;

    cout<<abs(-B)(50,25,5,2)<<endl;
    cout<<-B.abs()(50,25,5,2)<<endl;
    cout<<(-B).abs()(50,25,5,2)<<endl;
    cout<<endl<<endl<<endl;

    Array<dim3,dim4> D = (A[0][0]+B[0][0]*C[0][0]).abs() + 2;
    cout<<D<<endl;

    
    Array<dim1,dim2,dim3,dim4> K = ((A+B).eval() * C);
    cout<<K(50,25,5,2)<<endl;
    K = ((A+B) * C);
    cout<<K(50,25,5,2)<<endl;


    Array<dim1> M = ((A[0][4][5]+B[0][4][5]).eval() * C[0][4][5]);
    cout<<M(2)<<endl;
    M = ((A[0][4][5]+B[0][4][5]) * C[0][4][5]);
    cout<<M(2)<<endl;


    cout<<(A+B)[0][0]<<endl;
    cout<<(A+B).min()<<endl;
    cout<<min(A+B)<<endl;
    cout<<(A+B).max()<<endl;
    cout<<max(A+B)<<endl;


    Array<dim1> O(2);
    Array<dim1> P(3);
    cout<<(O+P)[0]<<endl;
    cout<<(O+P).min()<<endl;
    cout<<min(O+P)<<endl;
    cout<<(O+P).max()<<endl;
    cout<<max(O+P)<<endl;

    cout<<endl<<endl;

    cout<<(A==A)<<endl;
    cout<<(A==B)<<endl;
    cout<<(A!=A)<<endl;
    cout<<(A!=B)<<endl;

    cout<<endl<<endl<<endl;


    cout<<sin(A)(50,25,5,2)<<endl;
    cout<<cos(A)(50,25,5,2)<<endl;
    cout<<tan(A)(50,25,5,2)<<endl;
    cout<<exp(A)(50,25,5,2)<<endl;
    cout<<log(A)(50,25,5,2)<<endl;

    cout<<round(A+0.1)(50,25,5,2)<<endl;
    cout<<round(A+0.51)(50,25,5,2)<<endl;
    cout<<floor(A+0.1)(50,25,5,2)<<endl;
    cout<<ceil(A+0.1)(50,25,5,2)<<endl;




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