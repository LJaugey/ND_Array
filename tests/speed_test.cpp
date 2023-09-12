#include <valarray>
#include <chrono>
#include <iostream>

#include "../ND_Array/ND_Array.hpp"

typedef std::valarray<double> Mat_1;
typedef std::valarray<Mat_1> Mat_2;
typedef std::valarray<Mat_2> Mat_3;
typedef std::valarray<Mat_3> Mat_4;
typedef std::valarray<Mat_4> Mat_5;

using namespace std;


double t = 0;
const double dt = 0.001;
const double t_fin = 10;

const int N = 200;

const int dim = 3;

void func_par(ND::Array<double,dim,dim,dim> u)
{
    for(int l = 0; l<dim; l++)
    {
        for(int m = 0; m<dim; m++)
        {
            for(int s = 0; s<dim; s++)
            {
                u(l,m,s) += dt;
            }
        }
    }
}
void func_brac(ND::Array<double,dim,dim,dim> u)
{
    for(int l = 0; l<dim; l++)
    {
        for(int m = 0; m<dim; m++)
        {
            for(int s = 0; s<dim; s++)
            {
                u[l][m][s] += dt;
            }
        }
    }
}




void func_(Mat_3 & u)
{
    for(int l = 0; l<dim; l++)
    {
        for(int m = 0; m<dim; m++)
        {
            for(int s = 0; s<dim; s++)
            {
                u[l][m][s] += dt;
            }
        }
    }
}



int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    cout<<"================================"<<endl<<endl;
    cout<<"Testing lazy evaluation speed-up"<<endl<<endl;
    cout<<"================================"<<endl<<endl;

    ND::Array<double,N,N,dim,dim,dim> A(2.0);
    ND::Array<double,N,N,dim,dim,dim> B(0.5);
    ND::Array<double,N,N,dim,dim,dim> C(0.0);



    cout<<endl<<"Test start with lazy evaluation"<<endl;
    t = 0;
    start = std::chrono::high_resolution_clock::now();
    while(t<t_fin)
    {
        C = C + A*B*dt;
        t+=10*dt;
    }
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    cout<<endl<<"test end. Code ran in "<<(double)duration.count()/1000.0<<" seconds"<<endl<<endl<<endl;


    C = 0.0;
    
    cout<<endl<<"Test start with normal evaluation"<<endl;
    t = 0;
    start = std::chrono::high_resolution_clock::now();
    while(t<t_fin)
    {
        C = (C + (A*(B*dt).eval()).eval()).eval();
        t+=10*dt;
    }
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    cout<<endl<<"test end. Code ran in "<<(double)duration.count()/1000.0<<" seconds"<<endl<<endl<<endl;
    



    cout<<endl<<endl;


    ND::Array<double,N,N,dim,dim,dim> u;
    u.fill(1.0);

    Mat_5 u_(Mat_4(Mat_3(Mat_2(Mat_1(1.0, dim), dim), dim), N), N);


    cout<<"====================================="<<endl<<endl;
    cout<<"Testing operator() vs operator[] vs nested valarrays"<<endl<<endl;
    cout<<"====================================="<<endl<<endl;

    // Traverse the array once to reduce the difference with/without caching
    #pragma omp parallel for collapse(2)
    for(int i = 0; i<N; i++)
    {
        for(int j = 0; j<N; j++)
        {
            func_par(u[i][j]);
        }
    }

    cout<<endl<<"test start with operator()"<<endl;
    t = 0;
    start = std::chrono::high_resolution_clock::now();
    while(t<t_fin)
    {
        #pragma omp parallel for collapse(2)
        for(int i = 0; i<N; i++)
        {
            for(int j = 0; j<N; j++)
            {
                func_par(u[i][j]);
            }
        }
        t += dt;
    }
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    cout<<endl<<"test end. Code ran in "<<(double)duration.count()/1000.0<<" seconds"<<endl<<endl<<endl;
    
    
    cout<<endl<<"test start with operator[]"<<endl;
    t = 0;
    start = std::chrono::high_resolution_clock::now();
    while(t<t_fin)
    {
        #pragma omp parallel for collapse(2)
        for(int i = 0; i<N; i++)
        {
            for(int j = 0; j<N; j++)
            {
                func_brac(u[i][j]);
            }
        }
        t += dt;
    }
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    cout<<endl<<"test end. Code ran in "<<(double)duration.count()/1000.0<<" seconds"<<endl<<endl;

    
    t = 0;
    cout<<endl<<"test start with valarray"<<endl;
    start = std::chrono::high_resolution_clock::now();
    while(t<t_fin)
    {
        #pragma omp parallel for collapse(2)
        for(int i = 0; i<N; i++)
        {
            for(int j = 0; j<N; j++)
            {
                func_(u_[i][j]);
            }
        }
        t += dt;
    }
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    cout<<endl<<"test end. Code ran in "<<(double)duration.count()/1000.0<<" seconds"<<endl<<endl;

    return 0;
}