#ifndef LINEAR_ALGEMBRA_HPP
#define LINEAR_ALGEMBRA_HPP

#include "ND_array.hpp"


template<size_t dim1,size_t dim2>
class Matrix: public Array<dim1,dim2>
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

#endif