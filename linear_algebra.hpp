#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP

#include "ND_array.hpp"

#define TILE 32

namespace ND {

template<size_t dim1,size_t dim2>
class Matrix: public Array<dim1,dim2>
{
    using Array<dim1,dim2>::Array;
};



template<size_t dim>
class Vector: public Array<dim>
{
    using Array<dim>::Array;
};


template<size_t dim>
Matrix<dim,dim> diag(Vector<dim> V)
{
    Matrix<dim,dim> result;

    for(size_t i=0; i<dim; i++)
    {
        result(i,i) = V(i);
    }

    return result;
}


template<size_t dim1,size_t common_dim,size_t other_dim2>
Matrix<dim1,other_dim2> operator*(const Matrix<dim1,common_dim>& M1, const Matrix<common_dim,other_dim2>& M2)
{
    Matrix<dim1,other_dim2> result;
    
    if constexpr(dim1>TILE && common_dim>TILE)
    {
        constexpr size_t last_TILE_I = ((dim1-1)/TILE);
        constexpr size_t last_TILE_K = ((common_dim-1)/TILE);

        #pragma omp parallel for collapse(2) if(M1.length>PAR_SIZE)
        for(size_t I=0; I<last_TILE_I*TILE; I+=TILE)
        {
            for(size_t K=0; K<last_TILE_K*TILE; K+=TILE)
            {
                for(size_t i=I; i<I+TILE; i+=2)
                {
                    for(size_t k=K; k<K+TILE; k+=2)
                    {
                        for(size_t j=0; j<other_dim2; j++)
                        {
                            result(i+0,j) += M1(i+0,k+0)*M2(k+0,j) + M1(i+0,k+1)*M2(k+1,j);
                            result(i+1,j) += M1(i+1,k+0)*M2(k+0,j) + M1(i+1,k+1)*M2(k+1,j);
                        }
                    }
                }
            }
        }
        for(size_t I=0; I<dim1-TILE; I+=TILE)
        {
            for(size_t i=I; i<I+TILE; i+=2)
            {
                for(size_t k=last_TILE_K*TILE; k<common_dim; k++)
                {
                    for(size_t j=0; j<other_dim2; j++)
                    {
                        result(i+0,j) += M1(i+0,k+0)*M2(k+0,j);
                        result(i+1,j) += M1(i+1,k+0)*M2(k+0,j);
                    }
                }
            }
        }
        for(size_t K=0; K<common_dim-TILE; K+=TILE)
        {
            for(size_t i=last_TILE_I*TILE; i<dim1; i++)
            {
                for(size_t k=K; k<K+TILE; k+=2)
                {
                    for(size_t j=0; j<other_dim2; j++)
                    {
                        result(i+0,j) += M1(i+0,k+0)*M2(k+0,j) + M1(i+0,k+1)*M2(k+1,j);
                    }
                }
            }
        }
        for(size_t i=last_TILE_I*TILE; i<dim1; i++)
        {
            for(size_t k=last_TILE_K*TILE; k<common_dim; k++)
            {
                for(size_t j=0; j<other_dim2; j++)
                {
                    result(i+0,j) += M1(i+0,k+0)*M2(k+0,j);
                }
            }
        }
    }
    else if constexpr(dim1>TILE)
    {
        constexpr size_t last_TILE_I = ((dim1-1)/TILE);
        
        #pragma omp parallel for if(M1.length>PAR_SIZE)
        for(size_t I=0; I<last_TILE_I*TILE; I+=TILE)
        {
            for(size_t i=I; i<I+TILE; i+=2)
            {
                for(size_t k=0; k<common_dim; k++)
                {
                    for(size_t j=0; j<other_dim2; j++)
                    {
                        result(i+0,j) += M1(i+0,k+0)*M2(k+0,j);
                        result(i+1,j) += M1(i+1,k+0)*M2(k+0,j);
                    }
                }
            }
        }
        for(size_t i=last_TILE_I*TILE; i<dim1; i++)
        {
            for(size_t k=0; k<common_dim; k++)
            {
                for(size_t j=0; j<other_dim2; j++)
                {
                    result(i+0,j) += M1(i+0,k+0)*M2(k+0,j);
                }
            }
        }
    }
    else if constexpr(common_dim>TILE)
    {
        constexpr size_t last_TILE_K = ((common_dim-1)/TILE);

        #pragma omp parallel for if(M1.length>PAR_SIZE)
        for(size_t K=0; K<last_TILE_K*TILE; K+=TILE)
        {
            for(size_t i=0; i<dim1; i++)
            {
                for(size_t k=K; k<K+TILE; k+=2)
                {
                    for(size_t j=0; j<other_dim2; j++)
                    {
                        result(i+0,j) += M1(i+0,k+0)*M2(k+0,j) + M1(i+0,k+1)*M2(k+1,j);
                    }
                }
            }
        }
        for(size_t i=0; i<dim1; i++)
        {
            for(size_t k=last_TILE_K*TILE; k<common_dim; k++)
            {
                for(size_t j=0; j<other_dim2; j++)
                {
                    result(i+0,j) += M1(i+0,k+0)*M2(k+0,j);
                }
            }
        }
    }
    else
    {
        for(size_t i=0; i<dim1; i++)
        {
            for(size_t k=0; k<common_dim; k++)
            {
                for(size_t j=0; j<other_dim2; j++)
                {
                    result(i,j) += M1(i,k)*M2(k,j);
                }
            }
        }
    }
    
    return result;
}



template<size_t dim1,size_t common_dim>
Vector<dim1> operator*(const Matrix<dim1,common_dim>& M, const Vector<common_dim>& V)
{
    Vector<dim1> result;

    if constexpr(dim1>TILE)
    {
        constexpr size_t last_TILE_I = ((dim1-1)/TILE);
        
        #pragma omp parallel for if(M.length>PAR_SIZE)
        for(size_t I=0; I<last_TILE_I*TILE; I+=TILE)
        {
            for(size_t i=I; i<I+TILE; i+=2)
            {
                for(size_t j=0; j<common_dim; j++)
                {
                    result(i+0) += M(i+0,j)*V(j);
                    result(i+1) += M(i+1,j)*V(j);
                }
            }
        }
        for(size_t i=last_TILE_I*TILE; i<dim1; i++)
        {
            for(size_t j=0; j<common_dim; j++)
            {
                result(i+0) += M(i+0,j)*V(j);
            }
        }
    }
    else
    {
        #pragma omp parallel for if(M.length>PAR_SIZE)
        for(size_t i=0; i<dim1; i++)
        {
            for(size_t j=0; j<common_dim; j++)
            {
                result(i) += M(i,j)*V(j);
            }
        }
    }

    return result;
}


template<size_t common_dim>
double operator*(const Vector<common_dim>& V1, const Vector<common_dim>& V2)
{
    double result = 0.0;

    #pragma omp parallel for reduction(+:result) if(common_dim>PAR_SIZE)
    for(size_t i=0; i<common_dim; i++)
    {
        result += V1(i)*V2(i);
    }

    return result;
}

}

#endif
