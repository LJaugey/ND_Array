#ifndef ARRAY_EXPRESSION_HPP
#define ARRAY_EXPRESSION_HPP

#include <cstddef>

template <class E> 
struct base_traits;


template <class E>
class Array_Expression
{
public:
    
    typedef typename base_traits<E>::terminal_type terminal_type;

    
    inline double get_element(size_t i) const   {  return static_cast<const E&>(*this).get_element(i);   }

    template <typename... ind_type>
    inline double operator()(ind_type... indices)   {   return static_cast<const E&>(*this)(indices...);    }

    inline const E& abs() {  return abs(static_cast<const E&>(*this));   } 
};

#endif