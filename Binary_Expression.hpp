#ifndef BINARY_EXPRESSION_HPP
#define BINARY_EXPRESSION_HPP

#include <cstddef>
#include <type_traits>

#include "Array_Expression.hpp"
#include "Unary_Expression.hpp"


template <class E1, class OP, class E2>
class Binary_Op : public Array_Expression<Binary_Op<E1,OP,E2>>
{
    E1 arg1;
    E2 arg2;

public:
    
    typedef typename base_traits<Binary_Op>::terminal_type terminal_type;
    
    Binary_Op(E1 a_1, E2 a_2)
    :arg1(a_1),arg2(a_2)
    {}

    inline auto abs()   {   return Unary_Op<Array_abs,Binary_Op>(*this);    }

    inline auto get_element(size_t i) const
    {
        if constexpr(std::is_scalar_v<E1>)
        {
            return OP::apply(arg1,arg2.get_element(i));
        }
        else if constexpr(std::is_scalar_v<E2>)
        {
            return OP::apply(arg1.get_element(i),arg2);
        }
        else
        {
            return OP::apply(arg1.get_element(i),arg2.get_element(i));
        }
    }

    template <typename... ind_type>
    inline double operator()(ind_type... indices)
    {
        if constexpr(std::is_scalar_v<E1>)
        {
            return OP::apply(arg1,arg2(indices...));
        }
        else if constexpr(std::is_scalar_v<E2>)
        {
            return OP::apply(arg1(indices...),arg2);
        }
        else
        {
            return OP::apply(arg1(indices...),arg2(indices...));
        }
    }
};
template <class E1, class OP, class E2>
struct base_traits<Binary_Op<E1,OP,E2>>
{
    typedef typename std::conditional<  std::is_scalar_v<E1>,
                                        base_traits<E2>,
                                        base_traits<E1>
                                        >::type::terminal_type   terminal_type;
};


struct Array_add
{
    static inline double apply(double u,double v)   {   return u + v;  }
};
template <class LHS, class RHS>
auto operator+(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_add,RHS>(lhs,rhs);
}

struct Array_sub
{
    static inline double apply(double u, double v)  {   return u - v;  }
};
template <class LHS, class RHS>
auto operator-(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_sub,RHS>(lhs,rhs);
}

struct Array_mul
{
    static inline double apply(double u, double v)  {   return u * v;  }
};
template <class LHS, class RHS>
Binary_Op<LHS,Array_mul,RHS> operator*(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_mul,RHS>(lhs,rhs);
}

struct Array_div
{
    static inline double apply(double u, double v)  {   return u / v;  }
};
template <class LHS, class RHS>
Binary_Op<LHS,Array_div,RHS> operator/(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_div,RHS>(lhs,rhs);
}



#endif