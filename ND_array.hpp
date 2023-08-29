#ifndef NDARRAY_HPP
#define NDARRAY_HPP

#include <cstddef>
#include <iostream>
#include <type_traits>
#include <omp.h>


#include "Array_Expression.hpp"
#include "Unary_Expression.hpp"
#include "Binary_Expression.hpp"

#define PAR_SIZE 1024



template <size_t firstDim, size_t... RestDims>
class Array : public Array_Expression<Array<firstDim, RestDims...>>
{
public:

    static constexpr size_t N = sizeof...(RestDims) + 1;
    static constexpr size_t length = firstDim * (RestDims * ...);
    static constexpr size_t Dims[N] = {firstDim, RestDims...};

    typedef typename base_traits<Array>::terminal_type terminal_type;

protected:

    double* data_;
    bool is_original;


public:

    inline const double get_element(size_t i) const     {   return data_[i];    }
    inline double& get_element(size_t i)                {   return data_[i];    }

    // Base constructor
    Array()
    {
        data_ = new double[length];

        is_original =  true;
    }
    // Base constructor
    Array(double val)
    {
        data_ = new double[length];

        std::fill_n(data_,length, val);

        is_original =  true;
    }

    // copy constructor
    Array(const Array<firstDim, RestDims...>& other)
    {
        data_ = other.data_;
        is_original =  false;
    }

    // Constructor from pointer
    Array(double* p, bool is_or = false)
    {
        data_ = p;
        is_original = is_or;
    }
    
    // copy assigment operator
    const Array<firstDim, RestDims...>& operator=(const Array<firstDim, RestDims...>& other)
    {
        std::copy(other.data_, other.data_ + length, data_);
        return *this;
    }
    const Array<firstDim, RestDims...>& operator=(double val)
    {
        std::fill_n(data_,length, val);
        return *this;
    }
    // constructor from N-1 dimensional array
    Array(Array<RestDims...> const& slice)
    {
        data_ = new double[length];
        for(size_t i = 0; i<firstDim; i++)
        {
            std::copy(slice.data_, slice.data_ + slice.length, data_+ i*(RestDims*...));
        }
    }

    // construct from Array_expressions
    template <typename E>
    Array(Array_Expression<E> const& expr)
    : is_original(true)
    {
        data_ = new double[length];
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; ++i)
            {
                data_[i] = expr.get_element(i);
            }
        }
        else
        {
            for (size_t i = 0; i < length; ++i)
            {
                data_[i] = expr.get_element(i);
            }
        }
    }
    template <typename E>
    Array operator=(Array_Expression<E> const& expr)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; ++i)
            {
                data_[i] = expr.get_element(i);
            }
        }
        else
        {
            for (size_t i = 0; i < length; ++i)
            {
                data_[i] = expr.get_element(i);
            }
        }

        return *this;
    }

    // destructor
    ~Array() {  if(is_original)    delete[] data_;   }



    // access element
    template <typename... ind_type>
    inline double& operator()(ind_type... indices)
    {
        size_t offset = 0;
        size_t temp[N] = {static_cast<size_t>(indices)...};

        for (size_t i = 0; i < std::min(N, sizeof...(ind_type)); i++) {
            offset = offset * Dims[i] + temp[i];
        }
        return data_[offset];
    }
    template <typename... ind_type>
    inline const double& operator()(ind_type... indices) const
    {
        size_t offset = 0;
        size_t temp[N] = {static_cast<size_t>(indices)...};

        for (size_t i = 0; i < std::min(N, sizeof...(ind_type)); i++) {
            offset = offset * Dims[i] + temp[i];
        }
        return data_[offset];
    }


    // access element
    inline Array<RestDims...> operator[](size_t index)
    {
        return Array<RestDims...>(data_+ index * (RestDims * ...), false);
    }
    inline const Array<RestDims...> operator[](size_t index) const
    {
        return Array<RestDims...>(data_+ index * (RestDims * ...), false);
    }


    Array<firstDim, RestDims...>& fill(double val)
    {
        std::fill_n(data_,length, val);
        return *this;
    }
    Array<firstDim, RestDims...>& fill(const Array<firstDim, RestDims...>& other)
    {
        if(data_ != other.data_)    std::copy(other.data_, other.data_ + length, data_);
        return *this;
    }
    template<class E>
    Array<firstDim, RestDims...>& fill(Array_Expression<E> const& expr)
    {
        for (size_t i = 0; i < length; ++i)
        {
            data_[i] = expr[i];
        }

        return *this;
    }


    // explicit copy
    Array<firstDim, RestDims...> copy() const
    {
        Array<firstDim, RestDims...> result;

        std::copy(data_, data_+length, result.data_);

        return result;
    }


    inline const size_t size(const size_t index = 0) const
    {
        return Dims[index];
    }

    //comparison operators
    const bool operator==(Array<firstDim, RestDims...> const other) const
    {
        // No need to check the arrays' shape since this
        // method is only defined when the shape is the same
        
        for(size_t i=0; i<length; i++)
        {
            if(data_[i] != other.data_[i])  return false;
        }

        return true;
    }
    inline const bool operator!=(const Array<firstDim, RestDims...>& other) const   {   return !((*this)==other);   }

    template<class E>
    requires std::is_same_v<terminal_type, typename E::terminal_type>   //makes sure only arrays with the same shape are compared
    const bool operator==(const Array_Expression<E> expr) const
    {
        for(size_t i=0; i<length; i++)
        {
            if(data_[i] != expr.get_element(i)) return false;
        }

        return true;
    }
    template<class E>
    requires std::is_same_v<terminal_type, typename E::terminal_type>
    inline const bool operator!=(const Array_Expression<E> expr) const   {   return !((*this)==expr);   }




    // min
    const double min() const
    {
        double res = data_[0];
        
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for reduction(min:res) if(omp_get_num_threads() == 1)
            for (size_t i = 1; i < length; i++) {
                res = std::min(res,data_[i]);
            }
        }
        else
        {
            for (size_t i = 1; i < length; i++) {
                res = std::min(res,data_[i]);
            }
        }
        return res;
    }
    // max
    const double max() const
    {
        double res = data_[0];
        
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for reduction(min:res) if(omp_get_num_threads() == 1)
            for (size_t i = 1; i < length; i++) {
                res = std::max(res,data_[i]);
            }
        }
        else
        {
            for (size_t i = 1; i < length; i++) {
                res = std::max(res,data_[i]);
            }
        }
        return res;
    }





    // Arithmetic operations


    // += operator
    const Array<firstDim, RestDims...>& operator+=(const Array<firstDim, RestDims...>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] += other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] += other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const Array<firstDim, RestDims...>& operator+=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] += scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] += scalar;
            }
        }
        return *this;
    }

    // -= operator
    const Array<firstDim, RestDims...>& operator-=(const Array<firstDim, RestDims...>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] -= other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] -= other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const Array<firstDim, RestDims...>& operator-=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] -= scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] -= scalar;
            }
        }
        return *this;
    }

    // *= operator
    const Array<firstDim, RestDims...>& operator*=(const Array<firstDim, RestDims...>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const Array<firstDim, RestDims...>& operator*=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= scalar;
            }
        }
        return *this;
    }

    // /= operator
    const Array<firstDim, RestDims...>& operator/=(const Array<firstDim, RestDims...>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] /= other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] /= other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const Array<firstDim, RestDims...>& operator/=(double scalar)
    {
        double inv_scal = 1.0/scalar;

        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= inv_scal;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= inv_scal;
            }
        }
        return *this;
    }
};
template<size_t firstDim, size_t... RestDims>
struct base_traits<Array<firstDim, RestDims...>>
{
    typedef Array<firstDim, RestDims...> terminal_type;
};




// ostream
template <size_t firstDim, size_t... RestDims>
std::ostream& operator<<(std::ostream& output, const Array<firstDim, RestDims...>& other)
{
    if (sizeof...(RestDims) > 0)
    {
        for(size_t i = 0; i<firstDim; i++)
            output<<other[i]<<std::endl;
    }
    else
    {
        for(size_t i = 0; i<firstDim; i++)
            output<<other[i]<<"\t";
    }

    return output;
}



template <size_t Dim>
class Array<Dim>
{
public:

    static constexpr size_t N = 1;
    static constexpr size_t length = Dim;
    static constexpr size_t Dims[N] = {Dim};

    typedef typename base_traits<Array>::terminal_type terminal_type;

private:

    double* data_;
    bool is_original;

    
public:

    inline const double get_element(size_t i) const    {   return data_[i];    }
    inline double& get_element(size_t i)               {   return data_[i];    }

    // Base constructor
    Array()
    {
        data_ = new double[length];

        is_original =  true;
    }
    Array(double val)
    {
        data_ = new double[length];

        std::fill_n(data_,length, val);

        is_original =  true;
    }

    // copy constructor
    Array(const Array<Dim> & other)
    {
        data_ = other.data_;
        is_original =  false;
    }
    
    // Constructor from pointer
    Array(double* p, bool is_or = false)
    {
        data_ = p;
        is_original = is_or;
    }
    
    // copy assigment operator
    const Array<Dim>& operator=(const Array<Dim>& other)
    {
        std::copy(other.data_, other.data_ + length, data_);
        return *this;
    }
    const Array<Dim>& operator=(double val)
    {
        std::fill_n(data_,length, val);
        return *this;
    }
    
    // construct from Array_expressions
    template <typename E>
    Array(Array_Expression<E> const& expr)
    : is_original(true)
    {
        data_ = new double[length];
        for (size_t i = 0; i < length; ++i)
        {
            data_[i] = expr.get_element(i);
        }
    }
    template <typename E>
    Array operator=(Array_Expression<E> const& expr)
    {
        for (size_t i = 0; i < length; ++i)
        {
            data_[i] = expr.get_element(i);
        }

        return *this;
    }

    // destructor
    ~Array() {  if(is_original)    delete[] data_;   }



    // access element
    inline double& operator()(size_t index)                { return data_[index]; }
    inline const double operator()(size_t index) const     { return data_[index]; }
    inline double& operator[](size_t index)                { return data_[index]; }
    inline const double operator[](size_t index) const     { return data_[index]; }





    Array<Dim>& fill(double val)
    {
        std::fill_n(data_,length, val);
        return *this;
    }
    Array<Dim>& fill(const Array<Dim> & other)
    {
        if(data_ != other.data_)    std::copy(other.data_, other.data_ + length, data_);
        return *this;
    }
    template<class E>
    Array<Dim>& fill(Array_Expression<E> const& expr)
    {
        for (size_t i = 0; i < length; ++i)
        {
            data_[i] = expr[i];
        }

        return *this;
    }


    // explicit copy
    Array<Dim> copy() const
    {
        Array<Dim> result;

        std::copy(data_, data_+length, result.data_);

        return result;
    }





    // min
    const double min() const
    {
        double res = data_[0];
        
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for reduction(min:res) if(omp_get_num_threads() == 1)
            for (size_t i = 1; i < length; i++) {
                res = std::min(res,data_[i]);
            }
        }
        else
        {
            for (size_t i = 1; i < length; i++) {
                res = std::min(res,data_[i]);
            }
        }
        return res;
    }
    // max
    const double max() const
    {
        double res = data_[0];
        
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for reduction(min:res) if(omp_get_num_threads() == 1)
            for (size_t i = 1; i < length; i++) {
                res = std::max(res,data_[i]);
            }
        }
        else
        {
            for (size_t i = 1; i < length; i++) {
                res = std::max(res,data_[i]);
            }
        }
        return res;
    }





    // Arithmetic operations


    // += operator
    const Array<Dim>& operator+=(const Array<Dim>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] += other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] += other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const Array<Dim>& operator+=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] += scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] += scalar;
            }
        }
        return *this;
    }

    // -= operator
    const Array<Dim>& operator-=(const Array<Dim>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] -= other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] -= other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const Array<Dim>& operator-=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] -= scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] -= scalar;
            }
        }
        return *this;
    }

    // *= operator
    const Array<Dim>& operator*=(const Array<Dim>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const Array<Dim>& operator*=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= scalar;
            }
        }
        return *this;
    }

    // /= operator
    const Array<Dim>& operator/=(const Array<Dim>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] /= other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] /= other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const Array<Dim>& operator/=(double scalar)
    {
        double inv_scal = 1.0/scalar;

        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= inv_scal;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= inv_scal;
            }
        }
        return *this;
    }
};
template<size_t Dim>
struct base_traits<Array<Dim>>
{
    typedef Array<Dim> terminal_type;
};


#endif
