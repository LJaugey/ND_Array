#ifndef LINEAR_ALGEMBRA_HPP
#define LINEAR_ALGEMBRA_HPP

#include "ND_array.hpp"


template<size_t dim1,size_t dim2>
class Matrix: public Array<dim1,dim2>
{};

template<size_t dim>
class Vector: public Array<dim>
{};


template<size_t dim1,size_t common_dim,size_t other_dim2>
Matrix<dim1,other_dim2> operator*(const Matrix<dim1,common_dim>& M1, const Matrix<common_dim,other_dim2>& M2)
{
    Matrix<dim1,other_dim2> result;

    for(size_t i=0; i<dim1; i++)
    {
        for(size_t j=0; j<other_dim2; j++)
        {
            for(size_t k=0; k<common_dim; k++)
            {
                result(i,j) += M1(i,k)*M2(k,j);
            }
        }
    }

    return result;
}



template<size_t dim1,size_t common_dim>
Vector<dim1> operator*(const Matrix<dim1,common_dim>& M, const Vector<common_dim>& V)
{
    Vector<dim1> result;

    for(size_t i=0; i<dim1; i++)
    {
        for(size_t j=0; j<common_dim; j++)
        {
            result(i) += M(i,j)*V(j);
        }
    }

    return result;
}


template<size_t common_dim>
double operator*(const Vector<common_dim>& V1, const Vector<common_dim>& V2)
{
    double result = 0.0;

    for(size_t i=0; i<common_dim; i++)
    {
        result += V1(i)*V2(i);
    }

    return result;
}

#endif
