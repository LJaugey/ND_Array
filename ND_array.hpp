#ifndef NDARRAY_HPP
#define NDARRAY_HPP

#include <cstddef>
#include <iostream>
#include <type_traits>
#include <vector>
#include <omp.h>


#define PAR_SIZE 1024


template <class E>
class Array_Expression
{
public:
    
    inline double get_element(size_t i) const   {  return static_cast<const E&>(*this).get_element(i);   }

    inline size_t len() const   {   return static_cast<const E&>(*this).len();  }

    template <typename... ind_type>
    inline double operator()(ind_type... indices)              {   return static_cast<const E&>(*this)(indices...);    }
};

/*
// min
template<class E>
const double min(Array_Expression<E> Arr_expr)
{
    return Arr_expr.min();
}
// max
template<class E>
const double max(Array_Expression<E> Arr_expr)
{
    return Arr_expr.max();
}
*/


template <class E1, class OP, class E2>
class Binary_Op : public Array_Expression<Binary_Op<E1,OP,E2>>
{
    E1 arg1;
    E2 arg2;

    OP operation;

public:

    Binary_Op(E1 a_1, E2 a_2)
    :arg1(a_1),arg2(a_2)
    {}

    inline auto get_element(size_t i) const
    {
        if constexpr(std::is_convertible<E1, double>::value)
        {
            return operation.apply(arg1,arg2.get_element(i));
        }
        else if constexpr(std::is_convertible<E2, double>::value)
        {
            return operation.apply(arg1.get_element(i),arg2);
        }
        else
        {
            return operation.apply(arg1.get_element(i),arg2.get_element(i));
        }
    }

    inline size_t len() const
    {
        if constexpr(!std::is_convertible<E1, double>::value)
        {
            return arg1.len();
        }
        else if constexpr(!std::is_convertible<E2, double>::value)
        {
            return arg2.len();
        }
        else
        {
            return 1;
        }
    }

    template <typename... ind_type>
    inline double operator()(ind_type... indices)
    {
        if constexpr(std::is_convertible<E1, double>::value)
        {
            return operation.apply(arg1,arg2(indices...));
        }
        else if constexpr(std::is_convertible<E2, double>::value)
        {
            return operation.apply(arg1(indices...),arg2);
        }
        else
        {
            return operation.apply(arg1(indices...),arg2(indices...));
        }
    }
};

template <class OP, class E>
class Unary_Op : public Array_Expression<Unary_Op<OP,E>>
{
    E arg;

    OP operation;

public:

    Unary_Op(E a)
    :arg(a)
    {}

    inline auto get_element(size_t i) const
    {
        if constexpr(std::is_convertible<E, double>::value)
        {
            return operation.apply(arg);
        }
        else
        {
            return operation.apply(arg.get_element(i));
        }
    }

    inline size_t len() const
    {
        if constexpr(!std::is_convertible<E, double>::value)
        {
            return arg.len();
        }
        else
        {
            return 1;
        }
    }

    template <typename... ind_type>
    inline double operator()(ind_type... indices)
    {
        if constexpr(std::is_convertible<E, double>::value)
        {
            return operation.apply(arg);
        }
        else
        {
            return operation.apply(arg(indices...));
        }
    }
};



struct Array_add
{
    template<typename T1, typename T2>
    inline auto apply(T1 u,T2 v) const  {   return u + v;  }
};
template <class LHS, class RHS>
auto operator+(const LHS& lhs, const RHS& rhs) {
return Binary_Op<LHS,Array_add,RHS>(lhs,rhs);
}

struct Array_sub
{
    template<typename T1, typename T2>
    inline auto apply(T1 u,T2 v) const  {   return u - v;  }
};
template <class LHS, class RHS>
auto operator-(const LHS& lhs, const RHS& rhs) {
return Binary_Op<LHS,Array_sub,RHS>(lhs,rhs);
}

struct Array_mul
{
    template<typename T1, typename T2>
    inline auto apply(T1 u,T2 v) const  {   return u * v;  }
};
template <class LHS, class RHS>
Binary_Op<LHS,Array_mul,RHS> operator*(const LHS& lhs, const RHS& rhs) {
return Binary_Op<LHS,Array_mul,RHS>(lhs,rhs);
}

struct Array_div
{
    template<typename T1, typename T2>
    inline auto apply(T1 u,T2 v) const  {   return u / v;  }
};
template <class LHS, class RHS>
Binary_Op<LHS,Array_div,RHS> operator/(const LHS& lhs, const RHS& rhs) {
return Binary_Op<LHS,Array_div,RHS>(lhs,rhs);
}



struct Array_opp
{
    template<typename T>
    inline T apply(T u) const  {   return -u;  }
};
template <class RHS>
Unary_Op<Array_opp,RHS> operator-(const RHS& rhs) {
return Unary_Op<Array_opp,RHS>(rhs);
}

struct Array_abs
{
    template<typename T>
    inline T apply(T u) const  {   return std::abs(u);  }
};
template <class RHS>
Unary_Op<Array_opp,RHS> abs(const RHS& rhs) {
return Unary_Op<Array_opp,RHS>(rhs);
}







template <size_t firstDim, size_t... RestDims>
class Array : public Array_Expression<Array<firstDim, RestDims...>>
{
public:

    static constexpr size_t N = sizeof...(RestDims) + 1;
    static constexpr size_t length = firstDim * (RestDims * ...);
    static constexpr size_t Dims[N] = {firstDim, RestDims...};

protected:

    double* data_;
    bool is_original;


public:

    const double get_element(size_t i) const    {   return data_[i];    }
    double& get_element(size_t i)               {   return data_[i];    }
    inline size_t len() const                   {   return length;      }

    // Base constructor
    Array()
    {
        data_ = new double[length];
        
        std::fill_n(data_,length, 0.0);

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
            for (size_t i = 0; i < expr.len(); ++i)
            {
                data_[i] = expr.get_element(i);
            }
        }
        else
        {
            for (size_t i = 0; i < expr.len(); ++i)
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
            for (size_t i = 0; i < expr.len(); ++i)
            {
                data_[i] = expr.get_element(i);
            }
        }
        else
        {
            for (size_t i = 0; i < expr.len(); ++i)
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
        for (size_t i = 0; i < expr.len(); ++i)
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


    // abs
    Array<firstDim, RestDims...> const& abs()
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] = std::abs(data_[i]);
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] = std::abs(data_[i]);
            }
        }
        return *this;
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

private:

    double* data_;
    bool is_original;

    
public:

    const double get_element(size_t i) const    {   return data_[i];    }
    double& get_element(size_t i)               {   return data_[i];    }
    inline size_t len() const                   {   return length;      }

    // Base constructor
    Array()
    {
        data_ = new double[length];
        std::fill_n(data_,length, 0.0);

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
        for (size_t i = 0; i < expr.len(); ++i)
        {
            data_[i] = expr.get_element(i);
        }
    }
    template <typename E>
    Array operator=(Array_Expression<E> const& expr)
    {
        for (size_t i = 0; i < expr.len(); ++i)
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
        for (size_t i = 0; i < expr.len(); ++i)
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




    // abs
    Array<Dim> const& abs()
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] = std::abs(data_[i]);
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] = std::abs(data_[i]);
            }
        }
        return *this;
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


#endif
