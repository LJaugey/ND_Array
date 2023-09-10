#ifndef HELPER_HPP
#define HELPER_HPP

#include <stddef.h>
#include <type_traits>


#define PAR_SIZE 1024


namespace ND {


// Base traits
template <class E> 
struct base_traits;



// Array expression
template <class E>
class Array_Expression;


template <class E>
struct is_Array_Expression
{
    static constexpr bool value = std::is_base_of<Array_Expression<E>,E>::value;
};



// Unary operation
template <class OP, class E>
class Unary_Op;


template <class OP, class E>
struct base_traits<Unary_Op<OP,E>>
{
    typedef typename base_traits<E>::terminal_type terminal_type;

    typedef typename terminal_type::terminal_sub_type terminal_sub_type;

    typedef typename terminal_type::value_type value_type;
};



// Binary operation
template <class E1, class OP, class E2>
class Binary_Op;


template <class E1, class OP, class E2>
struct base_traits<Binary_Op<E1,OP,E2>>
{
    typedef typename std::conditional<  ND::is_Array_Expression<E1>::value,
                                        base_traits<E1>,
                                        base_traits<E2>
                                        >::type::terminal_type   terminal_type;

    typedef typename terminal_type::terminal_sub_type terminal_sub_type;
    
    typedef typename terminal_type::value_type value_type;
};



// ND Array
template<typename T, size_t firstDim, size_t... RestDims>
class Array;


template<typename T, size_t firstDim, size_t... RestDims>
struct base_traits<Array<T, firstDim, RestDims...>>
{
    typedef Array<T, firstDim, RestDims...> terminal_type;
    typedef Array<T, RestDims...> terminal_sub_type;
    typedef T value_type;
};


template<typename T, size_t Dim>
struct base_traits<Array<T, Dim>>
{
    typedef Array<T, Dim> terminal_type;
    typedef T terminal_sub_type;
    typedef T value_type;
};



}


#endif