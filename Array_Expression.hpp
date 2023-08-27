#ifndef ARRAY_EXPRESSION_HPP
#define ARRAY_EXPRESSION_HPP

#include <cstddef>
#include <iostream>
#include <type_traits>

template <class E> 
struct base_traits;

template<typename T>
requires std::is_scalar_v<T>
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

    inline const E& abs() {  return abs(static_cast<const E&>(*this));   } 


    terminal_type eval() const
    {
        return terminal_type(*this);
    }

    // For now, operator[] collapses the whole array
    const auto operator[](size_t index) const
    {
        return (this->eval())[index];
    }

    const double min() const
    {
        return (this->eval()).min();
    }
    const double max() const
    {
        return (this->eval()).max();
    }
};


template<class E>
const double min(const Array_Expression<E>& expr)
{
    return expr.min();
}
template<class E>
const double max(const Array_Expression<E>& expr)
{
    return expr.max();
}

template<class E>
std::ostream& operator<<(std::ostream& output, const Array_Expression<E>& expr)
{
    return output<<expr.eval();
}


#endif