#include <cstddef>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "../ND_Array/ND_Array.hpp"

const size_t dim_1 = 2;
const size_t dim_2 = 3;



TEST_CASE("Constructors") {

    const double a = 0.0;
    ND::Array<double, dim_2> A_sub(a);


    ND::Array<double,dim_1, dim_2> A1;
    ND::Array<double,dim_1, dim_2> A2(a);
    ND::Array<double,dim_1, dim_2> A3(0.0);
    ND::Array<double,dim_1, dim_2> A4(0);
    ND::Array<double,dim_1, dim_2> A5(A_sub);


    CHECK_EQ(A1,A2);
    CHECK_EQ(A1,A3);
    CHECK_EQ(A1,A4);
    CHECK_EQ(A1,A5);
    
    CHECK_EQ(A2,A3);
    CHECK_EQ(A2,A4);
    CHECK_EQ(A2,A5);
    
    CHECK_EQ(A3,A4);
    CHECK_EQ(A3,A5);
    
    CHECK_EQ(A4,A5);

    
    //copy constructor
    ND::Array<double,dim_1, dim_2> B(A1);
    ND::Array<double,dim_1, dim_2> B1;
    B1 = A1;

    CHECK_EQ(B,A1);
    CHECK_EQ(B1,A1);


    //from array expression
    const double b = 2.0;

    ND::Array<double,dim_1, dim_2> C(A1+b);
    ND::Array<double,dim_1, dim_2> C1;
    C1 = A1 + b;

    CHECK_EQ(C,(A1+b).eval());
    CHECK_EQ(C1,(A1+b).eval());
}


TEST_CASE("Comparison operations") {

    const double a = 1.0;

    ND::Array<double,dim_1, dim_2> A(a);
    ND::Array<double,dim_1, dim_2> A_2(2*a);

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

    ND::Array<double,dim_1, dim_2> A(a);
    ND::Array<double,dim_1, dim_2> B(b);
    ND::Array<double,dim_2> A_s(a);
    ND::Array<double,dim_2> B_s(b);

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

    ND::Array<double,dim_1, dim_2> A(a);
    ND::Array<double,dim_1, dim_2> B(b);
    ND::Array<double,dim_1, dim_2> C(c);



    SUBCASE("Addition") {
        ND::Array<double,dim_1, dim_2> A_B(a+b);
        ND::Array<double,dim_1, dim_2> A_B_C(a+b+c);

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
        ND::Array<double,dim_1, dim_2> A_B(a-b);
        ND::Array<double,dim_1, dim_2> A_B_C(a-b-c);

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
        ND::Array<double,dim_1, dim_2> A_B(a*b);
        ND::Array<double,dim_1, dim_2> A_B_C(a*b*c);

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
        ND::Array<double,dim_1, dim_2> A_B(a/b);
        ND::Array<double,dim_1, dim_2> A_B_C(a/b/c);

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

    
TEST_CASE("Trigonometric/Hyperbolic operations") {

    const double a = 0.4;
    const double b = -0.7;
    const double c = 1.7;

    ND::Array<double,dim_1, dim_2> A(a);
    ND::Array<double,dim_1, dim_2> B(b);
    ND::Array<double,dim_1, dim_2> C(c);


    
    SUBCASE("sin") {

        CHECK_EQ(sin(A)(0,0), sin(a));
        CHECK_EQ(sin(B)(0,0), sin(b));

        CHECK_EQ(asin(A)(0,0), asin(a));
        CHECK_EQ(asin(B)(0,0), asin(b));

        CHECK_EQ(sinh(A)(0,0), sinh(a));
        CHECK_EQ(sinh(B)(0,0), sinh(b));

        CHECK_EQ(asinh(A)(0,0), asinh(a));
        CHECK_EQ(asinh(B)(0,0), asinh(b));
    }
    
    SUBCASE("cos") {

        CHECK_EQ(cos(A)(0,0), cos(a));
        CHECK_EQ(cos(B)(0,0), cos(b));

        CHECK_EQ(acos(A)(0,0), acos(a));
        CHECK_EQ(acos(B)(0,0), acos(b));

        CHECK_EQ(cosh(A)(0,0), cosh(a));
        CHECK_EQ(cosh(B)(0,0), cosh(b));

        CHECK_EQ(acosh(C)(0,0), acosh(c));      // for acosh, arg>=1 is required
    }
    
    SUBCASE("tan") {

        CHECK_EQ(tan(A)(0,0), tan(a));
        CHECK_EQ(tan(B)(0,0), tan(b));

        CHECK_EQ(atan(A)(0,0), atan(a));
        CHECK_EQ(atan(B)(0,0), atan(b));

        CHECK_EQ(tanh(A)(0,0), tanh(a));
        CHECK_EQ(tanh(B)(0,0), tanh(b));

        CHECK_EQ(atanh(A)(0,0), atanh(a));
        CHECK_EQ(atanh(B)(0,0), atanh(b));
    }

    
    SUBCASE("tan2") {
        
        CHECK_EQ(atan2(A,B)(0,0), atan2(a,b));
        CHECK_EQ(atan2(B,A)(0,0), atan2(b,a));
        
        CHECK_EQ(atan2(A,b)(0,0), atan2(a,b));
        CHECK_EQ(atan2(a,B)(0,0), atan2(a,b));
    }
}



TEST_CASE("bitwise operations") {

    const int a = 5;
    const int b = 2;

    ND::Array<int,dim_1, dim_2> A(a);
    ND::Array<int,dim_1, dim_2> B(b);


    CHECK_EQ((A&B)(0,0), a&b);
    CHECK_EQ((A|B)(0,0), a|b);
    CHECK_EQ((A^B)(0,0), a^b);
    CHECK_EQ((A<<B)(0,0), a<<b);
    CHECK_EQ((A>>B)(0,0), a>>b);
}

    

TEST_CASE("Other") {

    const double a = 1.7;
    const double a_ = 1.2;

    ND::Array<double,dim_1, dim_2> A(a);
    ND::Array<double,dim_1, dim_2> A_(a_);

    CHECK_EQ(A.size(0), dim_1);
    CHECK_EQ(A.size(1), dim_2);


    CHECK_EQ((-A)(0,0), (-a));
    CHECK_EQ(abs(A)(0,0), std::abs(a));
    CHECK_EQ(exp(A)(0,0), std::exp(a));
    CHECK_EQ(log(A)(0,0), std::log(a));
    CHECK_EQ(log2(A)(0,0), std::log2(a));
    CHECK_EQ(log10(A)(0,0), std::log10(a));
    CHECK_EQ(sqrt(A)(0,0), std::sqrt(a));

    CHECK_EQ(pow(A,A_)(0,0), std::pow(a,a_));
    CHECK_EQ(pow(A,a_)(0,0), std::pow(a,a_));
    CHECK_EQ(pow(a_,A)(0,0), std::pow(a_,a));
    
    CHECK_EQ(round(A)(0,0), round(a));
    CHECK_EQ(ceil(A)(0,0), ceil(a));
    CHECK_EQ(floor(A)(0,0), floor(a));


    ND::Array<double,dim_1, dim_2> B;
    const double x = 1.0;
    const double y = 2.0;

    B[0] = x; B[1] = y;

    CHECK_EQ(B.min(), std::min(x,y));
    CHECK_EQ(min(B), std::min(x,y));

    CHECK_EQ(B.max(), std::max(x,y));
    CHECK_EQ(max(B), std::max(x,y));

    CHECK_EQ(B.sum(), (x+y)*dim_2);
    CHECK_EQ(sum(B), (x+y)*dim_2);


    int c = 5; int d = 2;
    ND::Array<int, dim_1,dim_2> C(c);
    ND::Array<int, dim_1,dim_2> D(d);

    CHECK_EQ((C%D)(0,0), c%d);
}