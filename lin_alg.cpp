#include <iostream>
#include "linear_algebra.hpp"
#include <chrono>


using namespace std;

#define SIZE 2500

int main()
{
    Matrix<SIZE,SIZE> A;
    A.fill(0.1);
    Matrix<SIZE,SIZE> B;
    B.fill(0.01);

    auto start = std::chrono::high_resolution_clock::now();
    Matrix<SIZE,SIZE> C = A*B;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    cout<<endl<<"Code ran in "<<(double)duration.count()/1000.0<<" seconds"<<endl<<endl;
    

    /*
    Matrix<12*SIZE,12*SIZE> M;

    Vector<12*SIZE> V;

    auto start = std::chrono::high_resolution_clock::now();
    Vector<12*SIZE> res = M*V;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    cout<<endl<<"Code ran in "<<(double)duration.count()/1000.0<<" seconds"<<endl<<endl;
    */

    /*
    Vector<SIZE*SIZE*80> A;
    A.fill(1.0/SIZE);
    Vector<SIZE*SIZE*80> B;
    B.fill(1.0/SIZE);

    auto start = std::chrono::high_resolution_clock::now();
    double res = A*B;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    cout<<endl<<"Code ran in "<<(double)duration.count()/1000.0<<" seconds"<<endl<<endl;
    */

    return 0;
}