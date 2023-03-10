/**
  * @file HermitePolynomials_test.cpp
  *
  * @brief 
  * Tests of functions to calculate probabilist Hermite polynomials 
  *
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "BasisFunctions.hpp" 
#include <gtest/gtest.h> 

TEST ( HermitePolynomial, H_0_ComplexInput ) {

    typedef std::complex<double> Complex;

    size_t index = 0;

    Complex x ( 3.0, 2.0 );
    Complex expected ( 1.0, 0.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_0_RealInput ) {

    size_t index = 0;

    double x = 7.0;
    double expected = 1.0;

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    EXPECT_FLOAT_EQ ( result, expected );

}

TEST ( HermitePolynomial, H_0_ImagInput ) {

    typedef std::complex<double> Complex;

    size_t index = 0;

    Complex x ( 0.0, 9.0 );
    Complex expected ( 1.0, 0.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_0_NegInput ) {

    typedef std::complex<double> Complex;

    size_t index = 0;

    Complex x ( -1.0, -2.0 );
    Complex expected ( 1.0, 0.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_0_ZeroInput ) {

    typedef std::complex<double> Complex;

    size_t index = 0;

    Complex x ( 0.0, 0.0 );
    Complex expected ( 1.0, 0.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_1_ComplexInput ) {

    typedef std::complex<double> Complex;

    size_t index = 1;

    Complex x ( 2.0, 2.0 );
    auto expected = x;

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_1_RealInput ) {

    size_t index = 1;

    double x = 8.0;
    double expected = x;

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    EXPECT_FLOAT_EQ ( result, expected );

}

TEST ( HermitePolynomial, H_1_ImagInput ) {

    typedef std::complex<double> Complex;

    size_t index = 1;

    Complex x ( 0.0, 2.0 );
    auto expected = x;

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_1_NegInput ) {

    typedef std::complex<double> Complex;

    size_t index = 1;

    Complex x ( 0.0, -5.0 );
    auto expected = x;

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_1_ZeroInput ) {

    typedef std::complex<double> Complex;

    size_t index = 1;

    Complex x ( 0.0, 0.0 );
    auto expected = x;

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_2_ComplexInput ) {

    typedef std::complex<double> Complex;

    size_t index = 2;

    Complex x ( 1.0, 3.0 );
    auto expected = x * x - 1.0;

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_3_ComplexInput ) {

    typedef std::complex<double> Complex;

    size_t index = 3;

    Complex x ( 4.0, 5.0 );
    auto expected = x * x * x - x * 3.0;

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_4_ComplexInput ) {

    typedef std::complex<double> Complex;

    size_t index = 4;

    Complex x ( -9.0, 2.0 );
    auto expected = x * x * x * x - x * x * 6.0 + 3.0;

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomials, H_0_0_0_0_4_ComplexInput ) {

    typedef std::complex<double> Complex;
    typedef std::vector<Complex> Vector;

    std::vector<size_t> indices { 0, 0, 0, 0 };

    Vector X {

        Complex (  2.0,  3.0 ), 
        Complex ( -5.0,  2.0 ), 
        Complex (  0.0,  2.0 ), 
        Complex ( -5.0,  0.0 ) 

    };

    size_t dim = 4;

    auto result =
        BasisFunctions::HermitePolynomials<

            size_t, double, std::complex<double>

        > ( indices, X, dim );

    Complex expected;
    expected.real(  1.0 );
    expected.imag(  0.0 );

    ASSERT_EQ ( result.size(), 1 );
    EXPECT_FLOAT_EQ ( result[0].real(), expected.real() );
    EXPECT_FLOAT_EQ ( result[0].imag(), expected.imag() );

}

TEST ( HermitePolynomials, H_2_0_0_0_4_ComplexInput ) {

    typedef std::complex<double> Complex;
    typedef std::vector<Complex> Vector;

    std::vector<size_t> indices { 2, 0, 0, 0 };

    Vector X {

        Complex (  2.0,  3.0 ), 
        Complex ( -5.0,  2.0 ), 
        Complex (  0.0,  2.0 ), 
        Complex ( -5.0,  0.0 ), 


    };

    size_t dim = 4;

    auto result =
        BasisFunctions::HermitePolynomials<

            size_t, double, std::complex<double>

        > ( indices, X, dim );

    Complex expected = ( X[0] * X[0] - 1.0 ) / std::sqrt(2);

    ASSERT_EQ ( result.size(), 1 );
    EXPECT_FLOAT_EQ ( result[0].real(), expected.real() );
    EXPECT_FLOAT_EQ ( result[0].imag(), expected.imag() );

}

TEST ( HermitePolynomials, H_1_0_1_0_4_ComplexInput ) {

    typedef std::complex<double> Complex;
    typedef std::vector<Complex> Vector;

    std::vector<size_t> indices { 1, 0, 1, 0 };

    Vector X {

        Complex (  2.0,  3.0 ), 
        Complex ( -5.0,  2.0 ), 
        Complex (  0.0,  2.0 ), 
        Complex ( -5.0,  0.0 ), 


    };

    size_t dim = 4;

    auto result =
        BasisFunctions::HermitePolynomials<

            size_t, double, std::complex<double>

        > ( indices, X, dim );

    Complex expected = X[0] * X[2];

    ASSERT_EQ ( result.size(), 1 );
    EXPECT_FLOAT_EQ ( result[0].real(), expected.real() );
    EXPECT_FLOAT_EQ ( result[0].imag(), expected.imag() );

}

TEST ( HermitePolynomials, H_0_0_0_3_4_ComplexInput ) {

    typedef std::complex<double> Complex;
    typedef std::vector<Complex> Vector;

    std::vector<size_t> indices { 0, 0, 0, 3 };

    Vector X {

        Complex (  2.0,  3.0 ), 
        Complex ( -5.0,  2.0 ), 
        Complex (  0.0,  2.0 ), 
        Complex ( -5.0,  0.0 ), 


    };

    size_t dim = 4;

    auto result =
        BasisFunctions::HermitePolynomials<

            size_t, double, std::complex<double>

        > ( indices, X, dim );

    Complex expected = ( X[3] * X[3] * X[3] - X[3] * 3.0 ) / std::sqrt(6);

    ASSERT_EQ ( result.size(), 1 );
    EXPECT_FLOAT_EQ ( result[0].real(), expected.real() );
    EXPECT_FLOAT_EQ ( result[0].imag(), expected.imag() );

}

TEST ( HermitePolynomials, H_2_1_1_2_0_3_2_ComplexInput ) {

    typedef std::complex<double> Complex;
    typedef std::vector<Complex> Vector;

    std::vector<size_t> indices { 2, 1, 1, 2, 0, 3 };

    Vector X {

        Complex (  2.0,  3.0 ), 
        Complex ( -5.0,  2.0 ) 

    };

    size_t dim = 2;

    auto result =
        BasisFunctions::HermitePolynomials<

            size_t, double, std::complex<double>

        > ( indices, X, dim );

    std::vector<std::complex<double>> expected;

    expected.push_back (
        X[1] * ( X[0] * X[0] - 1.0 ) / std::sqrt(2)
    );

    expected.push_back (
        X[0] * ( X[1] * X[1] - 1.0 ) / std::sqrt(2)
    );

    expected.push_back (
        ( X[1] * X[1] * X[1] - X[1] * 3.0 ) / std::sqrt(6)
    );

    ASSERT_EQ ( result.size(), expected.size() );

    for ( auto i = 0; i < expected.size(); i++ ) {
        
        EXPECT_FLOAT_EQ ( result[i].real(), expected[i].real() );
        EXPECT_FLOAT_EQ ( result[i].imag(), expected[i].imag() );

    }

}

TEST ( HermitePolynomials, H_Complex_ComplexInput ) {

    typedef std::complex<double> Complex;
    typedef std::vector<Complex> Vector;

    std::vector<size_t> indices { 2, 1, 1, 2, 0, 3 };

    Complex x1, x2, x3, x4;

    x1.real(  2.0 );
    x1.imag(  3.0 );

    x2.real( -5.0 );
    x2.imag(  2.0 );

    x3.real(  1.0 );
    x3.imag(  1.0 );

    x4.real( -1.0 );
    x4.imag(  2.0 );

    Vector X = { x1, x2, x3, x4 };

    size_t dim = 2;

    auto result =
        BasisFunctions::HermitePolynomials<

            size_t, double, std::complex<double>

        > ( indices, X, dim );

    Vector expected;

    expected.push_back (
        x2 * ( x1 * x1 - 1.0 ) / std::sqrt(2)
    );

    expected.push_back (
        x1 * ( x2 * x2 - 1.0 ) / std::sqrt(2)
    );

    expected.push_back (
        ( x2 * x2 * x2 - x2 * 3.0 ) / std::sqrt(6)
    );

    expected.push_back (
        x4 * ( x3 * x3 - 1.0 ) / std::sqrt(2)
    );

    expected.push_back (
        x3 * ( x4 * x4 - 1.0 ) / std::sqrt(2)
    );

    expected.push_back (
        ( x4 * x4 * x4 - x4 * 3.0 ) / std::sqrt(6)
    );

    ASSERT_EQ ( result.size(), expected.size() );

    for ( auto i = 0; i < expected.size(); i++ ) {
        
        EXPECT_FLOAT_EQ ( result[i].real(), expected[i].real() );
        EXPECT_FLOAT_EQ ( result[i].imag(), expected[i].imag() );

    }

}


TEST ( HermitePolynomials, WrongXSize ) {

    typedef std::complex<double> Complex;
    typedef std::vector<Complex> Vector;

    std::vector<size_t> indices { 2, 1, 1, 2, 0, 3 };

    Complex x1, x2;

    x1.real(  2.0 );
    x1.imag(  3.0 );

    x2.real( -5.0 );
    x2.imag(  2.0 );

     Vector X = { x1, x2 };

    size_t dim = 3;

    try {

        auto result =
            BasisFunctions::HermitePolynomials<

                size_t, double, std::complex<double>

            > ( indices, X, dim );

    } catch ( std::exception& e ) {

        EXPECT_STREQ (
            "HermitePolynomials: num of samples not multiple of dimension",
            e.what()
        );

    }

}

TEST ( HermitePolynomials, WrongIndicesSize ) {

    typedef std::complex<double> Complex;
    typedef std::vector<Complex> Vector;

    std::vector<size_t> indices { 2, 1, 1, 2, 0 };

    Complex x1, x2;

    x1.real(  2.0 );
    x1.imag(  3.0 );

    x2.real( -5.0 );
    x2.imag(  2.0 );

    Vector X = { x1, x2 };

    size_t dim = 2;

    try {

        auto result =
            BasisFunctions::HermitePolynomials<

                size_t, double, std::complex<double>

            > ( indices, X, dim );

    } catch ( std::exception& e ) {

        EXPECT_STREQ (
            "HermitePolynomials: num of indices not multiple of dimension", 
            e.what()
        );

    }

}

TEST ( HermitePolynomials, ZeroDimension ) {

    typedef std::complex<double> Complex;
    typedef std::vector<Complex> Vector;

    std::vector<size_t> indices { 2, 1, 1, 2, 0, 3 };

    Complex x1, x2;

    x1.real(  2.0 );
    x1.imag(  3.0 );

    x2.real( -5.0 );
    x2.imag(  2.0 );

    Vector X = { x1, x2 };

    size_t dim = 0;

    try {

        auto result =
            BasisFunctions::HermitePolynomials<

                size_t, double, std::complex<double>

            > ( indices, X, dim );

    } catch ( std::exception& e ) {

        EXPECT_STREQ (
            "HermitePolynomials: dimension must be positive", 
            e.what()
        );

    }

}

