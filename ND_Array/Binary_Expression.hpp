#ifndef BINARY_EXPRESSION_HPP
#define BINARY_EXPRESSION_HPP


#include "helper.hpp"
#include "Array_Expression.hpp"
#include "Unary_Expression.hpp"

namespace ND {

template <class E1, class OP, class E2>
class Binary_Op : public Array_Expression<Binary_Op<E1,OP,E2>>
{
    const E1& arg1;
    const E2& arg2;

public:
    
    typedef typename base_traits<Binary_Op>::terminal_type terminal_type;
    typedef typename base_traits<Binary_Op>::terminal_sub_type terminal_sub_type;
    typedef typename base_traits<Binary_Op>::value_type value_type;
    
    Binary_Op(const E1& a_1, const E2& a_2)
    :arg1(a_1),arg2(a_2)
    {}

    inline const value_type get_element(const size_t i) const
    {
        if constexpr(not ND::is_Array_Expression<E1>::value)
        {
            return OP::apply(arg1,arg2.get_element(i));
        }
        else if constexpr(not ND::is_Array_Expression<E2>::value)
        {
            return OP::apply(arg1.get_element(i),arg2);
        }
        else
        {
            return OP::apply(arg1.get_element(i),arg2.get_element(i));
        }
    }

    template <typename... ind_type>
    inline const value_type operator()(const ind_type... indices) const
    {
        if constexpr(not ND::is_Array_Expression<E1>::value)
        {
            return OP::apply(arg1,arg2(indices...));
        }
        else if constexpr(not ND::is_Array_Expression<E2>::value)
        {
            return OP::apply(arg1(indices...),arg2);
        }
        else
        {
            return OP::apply(arg1(indices...),arg2(indices...));
        }
    }
};




struct Array_add
{
    static inline const auto apply(const auto u, const auto v)   {   return u + v;  }
};
template <class LHS, class RHS>
static inline const Binary_Op<LHS,Array_add,RHS> operator+(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_add,RHS>(lhs,rhs);
}

struct Array_sub
{
    static inline const auto apply(const auto u, const auto v)  {   return u - v;  }
};
template <class LHS, class RHS>
static inline const Binary_Op<LHS,Array_sub,RHS> operator-(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_sub,RHS>(lhs,rhs);
}

struct Array_mul
{
    static inline const auto apply(const auto u, const auto v)  {   return u * v;  }
};
template <class LHS, class RHS>
static inline const Binary_Op<LHS,Array_mul,RHS> operator*(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_mul,RHS>(lhs,rhs);
}

struct Array_div
{
    static inline const auto apply(const auto u, const auto v)  {   return u / v;  }
};
template <class LHS, class RHS>
static inline const Binary_Op<LHS,Array_div,RHS> operator/(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_div,RHS>(lhs,rhs);
}

struct Array_mod
{
    static inline const auto apply(const auto u, const auto v)  {   return u % v;  }
};
template <class LHS, class RHS>
static inline const Binary_Op<LHS,Array_mod,RHS> operator%(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_mod,RHS>(lhs,rhs);
}

struct Array_and
{
    static inline const auto apply(const auto u, const auto v)  {   return u & v;  }
};
template <class LHS, class RHS>
static inline const Binary_Op<LHS,Array_and,RHS> operator&(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_and,RHS>(lhs,rhs);
}

struct Array_or
{
    static inline const auto apply(const auto u, const auto v)  {   return u | v;  }
};
template <class LHS, class RHS>
static inline const Binary_Op<LHS,Array_or,RHS> operator|(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_or,RHS>(lhs,rhs);
}

struct Array_xor
{
    static inline const auto apply(const auto u, const auto v)  {   return u ^ v;  }
};
template <class LHS, class RHS>
static inline const Binary_Op<LHS,Array_xor,RHS> operator^(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_xor,RHS>(lhs,rhs);
}

struct Array_lshift
{
    static inline const auto apply(const auto u, const auto v)  {   return u << v;  }
};
template <class LHS, class RHS>
static inline const Binary_Op<LHS,Array_lshift,RHS> operator<<(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_lshift,RHS>(lhs,rhs);
}

struct Array_rshift
{
    static inline const auto apply(const auto u, const auto v)  {   return u >> v;  }
};
template <class LHS, class RHS>
static inline const Binary_Op<LHS,Array_rshift,RHS> operator>>(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_rshift,RHS>(lhs,rhs);
}



struct Array_pow
{
    static inline const auto apply(const auto u, const auto v)  {   return std::pow(u,v);  }
};
template <class LHS, class RHS>
static inline const Binary_Op<LHS,Array_pow,RHS> operator||(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_pow,RHS>(lhs,rhs);
}



struct Array_atan2
{
    static inline const auto apply(const auto u, const auto v)  {   return atan2(u,v);  }
};
template <class LHS, class RHS>
static inline const Binary_Op<LHS,Array_atan2,RHS> atan2(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_atan2,RHS>(lhs,rhs);
}



}

#endif