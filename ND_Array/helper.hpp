#ifndef HELPER_HPP
#define HELPER_HPP

#include <stddef.h>
#include <type_traits>


#define PAR_SIZE 1024


#ifdef _OPENMP
    #include <omp.h>
    #define STRINGIFY(a) #a
    #define OMP_FOR(n) _Pragma(STRINGIFY(omp parallel for if(n>PAR_SIZE)))
    #define OMP_FOR_add(n,var) _Pragma(STRINGIFY(omp parallel for reduction(+:var) if(n>PAR_SIZE)))
    #define OMP_FOR_min(n,var) _Pragma(STRINGIFY(omp parallel for reduction(min:var) if(n>PAR_SIZE)))
    #define OMP_FOR_max(n,var) _Pragma(STRINGIFY(omp parallel for reduction(max:var) if(n>PAR_SIZE)))
#else
    #define omp_get_thread_num() 0
    #define OMP_FOR(n)
    #define OMP_FOR_add(n,var)
    #define OMP_FOR_min(n,var)
    #define OMP_FOR_max(n,var)
#endif



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