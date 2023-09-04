#include <cstddef>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "../ND_array.hpp"


const size_t dim_1 = 2;
const size_t dim_2 = 3;


TEST_CASE("Comparison operations") {

    const double a = 1.0;

    ND::Array<dim_1, dim_2> A(a);
    ND::Array<dim_1, dim_2> A_2(2*a);

    CHECK_EQ(A,A);
    CHECK_EQ(A_2,A+A);
    CHECK_EQ(A+A,A_2);
    CHECK_EQ(A+A,A+A);

    CHECK(A!=A+A);
    CHECK(A+A!=A);
    CHECK_FALSE(A==A+A);
}

TEST_CASE("Access operators") {

    const double a = 2.0;
    const double b = 1.0;

    ND::Array<dim_1, dim_2> A(a);
    ND::Array<dim_1, dim_2> B(b);
    ND::Array<dim_2> A_s(a);
    ND::Array<dim_2> B_s(b);

    CHECK_EQ(A,A);
    CHECK_EQ(A[0],A_s);
    CHECK_EQ(A[dim_1-1],A_s);
    CHECK_EQ(A[0][0],a);

    CHECK_EQ(A(0,0),a);
    CHECK_EQ(A(dim_1-1,dim_2-1),a);

    A(0,0) = b;
    CHECK_EQ(A(0,0),b);
    CHECK_EQ(A(dim_1-1,dim_2-1),a);
    A(0,0) = a;

    CHECK_EQ((A+B)(0,0),a+b);
    CHECK_EQ((A+b)(0,0),a+b);
    CHECK_EQ((A+b)(dim_1-1,dim_2-1),a+b);
    CHECK_EQ((A+B)[0],A_s+B_s);
    CHECK_EQ((A+b)[0][0],a+b);
    CHECK_EQ((A+b)[dim_1-1][dim_2-1],a+b);
    
    A[0] = b;
    CHECK_EQ(A[0],B_s);
    CHECK_FALSE(A[dim_1-1]==B_s);
    CHECK_EQ(A[dim_1-1],A_s);
    CHECK_EQ(A[0][0],b);

    CHECK_EQ(A(0,0),b);
    CHECK_EQ(A(dim_1-1,dim_2-1),a);
}

TEST_CASE("Arithmetic operations") {

    const double a = 1.0;
    const double b = 0.2;
    const double c = -0.7;

    ND::Array<dim_1, dim_2> A(a);
    ND::Array<dim_1, dim_2> B(b);
    ND::Array<dim_1, dim_2> C(c);



    SUBCASE("Addition") {
        ND::Array<dim_1, dim_2> A_B(a+b);
        ND::Array<dim_1, dim_2> A_B_C(a+b+c);

        CHECK_EQ(A+B, A_B);
        CHECK_EQ(A+b, A_B);
        CHECK_EQ(a+B, A_B);

        CHECK_EQ(A+B+C, A_B_C);
        CHECK_EQ(A+B+c, A_B_C);
        CHECK_EQ(A+b+C, A_B_C);
        CHECK_EQ(a+B+C, A_B_C);

        CHECK_EQ(A+=B+C, a+(B+C));
        CHECK_EQ(B+=c, b+C);

    }
    SUBCASE("Subtraction") {
        ND::Array<dim_1, dim_2> A_B(a-b);
        ND::Array<dim_1, dim_2> A_B_C(a-b-c);

        CHECK_EQ(A-B, A_B);
        CHECK_EQ(A-b, A_B);
        CHECK_EQ(a-B, A_B);

        CHECK_EQ(A-B-C, A_B_C);
        CHECK_EQ(A-B-c, A_B_C);
        CHECK_EQ(A-b-C, A_B_C);
        CHECK_EQ(a-B-C, A_B_C);

        CHECK_EQ(A-=B-C, a-(B-C));
        CHECK_EQ(B-=c, b-C);
    }
    SUBCASE("Multiplication") {
        ND::Array<dim_1, dim_2> A_B(a*b);
        ND::Array<dim_1, dim_2> A_B_C(a*b*c);

        CHECK_EQ(A*B, A_B);
        CHECK_EQ(A*b, A_B);
        CHECK_EQ(a*B, A_B);

        CHECK_EQ(A*B*C, A_B_C);
        CHECK_EQ(A*B*c, A_B_C);
        CHECK_EQ(A*b*C, A_B_C);
        CHECK_EQ(a*B*C, A_B_C);

        CHECK_EQ(A*=B*C, a*(B*C));
        CHECK_EQ(B*=c, b*C);
    }
    SUBCASE("Division") {
        ND::Array<dim_1, dim_2> A_B(a/b);
        ND::Array<dim_1, dim_2> A_B_C(a/b/c);

        CHECK_EQ(A/B, A_B);
        CHECK_EQ(A/b, A_B);
        CHECK_EQ(a/B, A_B);

        CHECK_EQ(A/B/C, A_B_C);
        CHECK_EQ(A/B/c, A_B_C);
        CHECK_EQ(A/b/C, A_B_C);
        CHECK_EQ(a/B/C, A_B_C);

        CHECK_EQ(A/=B/C, a/(B/C));
        CHECK_EQ(B/=c, b/C);
    }
}

