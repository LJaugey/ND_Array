#ifndef ARRAY_EXPRESSION_HPP
#define ARRAY_EXPRESSION_HPP

#include <cstddef>
#include <iostream>
#include <type_traits>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#define PAR_SIZE 1024



namespace ND {

template <class E> 
struct base_traits;

template<typename T>
requires std::is_arithmetic_v<T>
struct base_traits<T>
{
    typedef double terminal_type;
};


template <class E>
class Array_Expression
{
public:
    
    typedef typename base_traits<E>::terminal_type terminal_type;

    
    inline const double get_element(size_t i) const   {  return static_cast<const E&>(*this).get_element(i);   }

    template <typename... ind_type>
    inline const double operator()(ind_type... indices) const   {   return static_cast<const E&>(*this)(indices...);    }


    inline const terminal_type eval() const
    {
        return terminal_type(*this);    // Guaranteed copy elision
    }

    // For simplicity, operator[] collapses the whole array
    // Could compute only the sub-array arr[i].
    // Could be a problem when N=1
    inline const auto operator[](size_t index) const
    {
        return (this->eval())[index];
    }

    



    // min
    const double min() const
    {
        double res = get_element(0);
        
        if constexpr(terminal_type::length>PAR_SIZE)
        {
            #pragma omp parallel for reduction(min:res) if(omp_get_num_threads() == 1)
            for (size_t i = 1; i < terminal_type::length; i++)
                res = std::min(res,get_element(i));
        }
        else {
            for (size_t i = 1; i < terminal_type::length; i++)
                res = std::min(res,get_element(i));
        }
        return res;
    }
    // max
    const double max() const
    {
        double res = get_element(0);
        
        if constexpr(terminal_type::length>PAR_SIZE)
        {
            #pragma omp parallel for reduction(max:res) if(omp_get_num_threads() == 1)
            for (size_t i = 1; i < terminal_type::terminal_type::length; i++)
                res = std::max(res,get_element(i));
        }
        else {
            for (size_t i = 1; i < terminal_type::terminal_type::length; i++) {
                res = std::max(res,get_element(i));
            }
        }
        return res;
    }
    // sum
    const double sum() const
    {
        double res = get_element(0);
        
        if constexpr(terminal_type::terminal_type::length>PAR_SIZE)
        {
            #pragma omp parallel for reduction(+:res) if(omp_get_num_threads() == 1)
            for (size_t i = 1; i < terminal_type::terminal_type::length; i++)
                res += get_element(i);
        }
        else {
            for (size_t i = 1; i < terminal_type::terminal_type::length; i++) {
                res += get_element(i);
            }
        }
        return res;
    }
    // mean
    inline const double mean() const   {   return (this->sum())/terminal_type::length; }

    const double stdev() const
    {
        double res = 0.0;
        double m = this->mean();
        
        if constexpr(terminal_type::length>PAR_SIZE)
        {
            #pragma omp parallel for reduction(+:res) if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < terminal_type::length; i++)
                res += get_element(i)*get_element(i);
        }
        else {
            for (size_t i = 0; i < terminal_type::length; i++)
                res += get_element(i)*get_element(i);
        }
        return sqrt(res/terminal_type::length - m*m);
    }

};


template<class E>
inline const double min(const Array_Expression<E>& expr)   {   return expr.min();  }
template<class E>
inline const double max(const Array_Expression<E>& expr)   {   return expr.max();  }
template<class E>
inline const double sum(const Array_Expression<E>& expr)   {   return expr.sum();  }
template<class E>
inline const double mean(const Array_Expression<E>& expr)  {   return expr.mean(); }
template<class E>
inline const double stdev(const Array_Expression<E>& expr) {   return expr.stdev();}

template<class E>
inline std::ostream& operator<<(std::ostream& output, const Array_Expression<E>& expr)
{
    return output<<expr.eval();
}


//comparison operators
template<class E1, class E2>
const bool operator==(const Array_Expression<E1>& expr1, const Array_Expression<E2>& expr2)
{
    if constexpr(!std::is_same_v<typename E1::terminal_type, typename E2::terminal_type>)  return false;

    for(size_t i=0; i<E1::terminal_type::length; i++)
    {
        if(expr1.get_element(i) != expr2.get_element(i))  return false;
    }

    return true;
}
template<class E1, class E2>
inline const bool operator!=(const Array_Expression<E1>& expr1, const Array_Expression<E2>& expr2)
{
    return !(expr1==expr2);
}


}

#endif