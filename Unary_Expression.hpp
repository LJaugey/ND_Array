#ifndef UNARY_EXPRESSION_HPP
#define UNARY_EXPRESSION_HPP

#include <cstddef>
#include <type_traits>
#include <algorithm>

#include "Array_Expression.hpp"


template <class OP, class E>
class Unary_Op : public Array_Expression<Unary_Op<OP,E>>
{
    E arg;

public:

    typedef typename base_traits<E>::terminal_type terminal_type;


    Unary_Op(E a)
    :arg(a)
    {}

    inline auto abs()
    {
        return arg.abs();
    }

    inline auto get_element(size_t i) const
    {
        if constexpr(std::is_convertible<E, double>::value)
        {
            return OP::apply(arg);
        }
        else
        {
            return OP::apply(arg.get_element(i));
        }
    }

    template <typename... ind_type>
    inline double operator()(ind_type... indices)
    {
        if constexpr(std::is_convertible<E, double>::value)
        {
            return OP::apply(arg);
        }
        else
        {
            return OP::apply(arg(indices...));
        }
    }
};
template <class OP, class E>
struct base_traits<Unary_Op<OP,E>>
{
    typedef typename base_traits<E>::terminal_type terminal_type;
};



struct Array_opp
{
    static inline double apply(double u)    {   return -u;  }
};
template <class RHS>
Unary_Op<Array_opp,RHS> operator-(const RHS& rhs)
{
    return Unary_Op<Array_opp,RHS>(rhs);
}

struct Array_abs
{
    static inline double apply(double u)    {   return std::abs(u);  }
};
template <class RHS>
Unary_Op<Array_abs,RHS> abs(const RHS& rhs)
{
    return Unary_Op<Array_abs,RHS>(rhs);
}



#endif