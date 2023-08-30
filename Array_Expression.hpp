#ifndef ARRAY_EXPRESSION_HPP
#define ARRAY_EXPRESSION_HPP

#include <cstddef>
#include <iostream>
#include <type_traits>

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

    
    inline double get_element(size_t i) const   {  return static_cast<const E&>(*this).get_element(i);   }

    template <typename... ind_type>
    inline double operator()(ind_type... indices)   {   return static_cast<const E&>(*this)(indices...);    }


    const terminal_type eval() const
    {
        return terminal_type(*this);
    }

    // For simplicity, operator[] collapses the whole array
    // Could compute only the sub-array arr[i].
    // Could be a problem when N=1
    const auto operator[](size_t index) const
    {
        return (this->eval())[index];
    }

    const double min() const    {   return (this->eval()).min();    }
    const double max() const    {   return (this->eval()).max();    }
    const double sum() const    {   return (this->eval()).sum();    }
    const double mean() const   {   return (this->eval()).mean();   }
    const double stdev() const  {   return (this->eval()).stdev();  }
};


template<class E>
const double min(const Array_Expression<E>& expr)   {   return expr.min();  }
template<class E>
const double max(const Array_Expression<E>& expr)   {   return expr.max();  }
template<class E>
const double sum(const Array_Expression<E>& expr)   {   return expr.sum();  }
template<class E>
const double mean(const Array_Expression<E>& expr)  {   return expr.mean(); }
template<class E>
const double stdev(const Array_Expression<E>& expr) {   return expr.stdev();}

template<class E>
std::ostream& operator<<(std::ostream& output, const Array_Expression<E>& expr)
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