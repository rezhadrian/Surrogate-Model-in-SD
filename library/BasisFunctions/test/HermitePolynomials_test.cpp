/**
  * @file HermitePolynomials_test.cpp
  *
  * @brief test the function HermitePolynomials 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "BasisFunctions.hpp" 
#include <gtest/gtest.h> 

TEST ( HermitePolynomials, H_0_0_0_0_4_ComplexInput ) {

    std::vector<size_t> indices { 0, 0, 0, 0 };

    std::complex<double> x1, x2, x3, x4;

    x1.real(  2.0 );
    x1.imag(  3.0 );

    x2.real( -5.0 );
    x2.imag(  2.0 );
        
    x3.real(  0.0 );
    x3.imag(  2.0 );

    x4.real( -5.0 );
    x4.imag(  0.0 );

    std::vector<std::complex<double>> X = { x1, x2, x3, x4 };

    size_t dim = 4;

    auto result =
        BasisFunctions::HermitePolynomials ( indices, X, dim );

    std::complex<double> expected;
    expected.real(  1.0 );
    expected.imag(  0.0 );

    ASSERT_EQ ( result.size(), 1 );
    EXPECT_FLOAT_EQ ( result[0].real(), expected.real() );
    EXPECT_FLOAT_EQ ( result[0].imag(), expected.imag() );

}

TEST ( HermitePolynomials, H_2_0_0_0_4_ComplexInput ) {

    std::vector<size_t> indices { 2, 0, 0, 0 };

    std::complex<double> x1, x2, x3, x4;

    x1.real(  2.0 );
    x1.imag(  3.0 );

    x2.real( -5.0 );
    x2.imag(  2.0 );
        
    x3.real(  0.0 );
    x3.imag(  2.0 );

    x4.real( -5.0 );
    x4.imag(  0.0 );

    std::vector<std::complex<double>> X = { x1, x2, x3, x4 };

    size_t dim = 4;

    auto result =
        BasisFunctions::HermitePolynomials ( indices, X, dim );

    std::complex<double> expected = ( x1 * x1 - 1.0 ) / std::sqrt(2);

    ASSERT_EQ ( result.size(), 1 );
    EXPECT_FLOAT_EQ ( result[0].real(), expected.real() );
    EXPECT_FLOAT_EQ ( result[0].imag(), expected.imag() );

}

TEST ( HermitePolynomials, H_1_0_1_0_4_ComplexInput ) {

    std::vector<size_t> indices { 1, 0, 1, 0 };

    std::complex<double> x1, x2, x3, x4;

    x1.real(  2.0 );
    x1.imag(  3.0 );

    x2.real( -5.0 );
    x2.imag(  2.0 );
        
    x3.real(  0.0 );
    x3.imag(  2.0 );

    x4.real( -5.0 );
    x4.imag(  0.0 );

    std::vector<std::complex<double>> X = { x1, x2, x3, x4 };

    size_t dim = 4;

    auto result =
        BasisFunctions::HermitePolynomials ( indices, X, dim );

    std::complex<double> expected = x1 * x3;

    ASSERT_EQ ( result.size(), 1 );
    EXPECT_FLOAT_EQ ( result[0].real(), expected.real() );
    EXPECT_FLOAT_EQ ( result[0].imag(), expected.imag() );

}

TEST ( HermitePolynomials, H_0_0_0_3_4_ComplexInput ) {

    std::vector<size_t> indices { 0, 0, 0, 3 };

    std::complex<double> x1, x2, x3, x4;

    x1.real(  2.0 );
    x1.imag(  3.0 );

    x2.real( -5.0 );
    x2.imag(  2.0 );
        
    x3.real(  0.0 );
    x3.imag(  2.0 );

    x4.real( -5.0 );
    x4.imag(  0.0 );

    std::vector<std::complex<double>> X = { x1, x2, x3, x4 };

    size_t dim = 4;

    auto result =
        BasisFunctions::HermitePolynomials ( indices, X, dim );

    std::complex<double> expected = ( x4 * x4 * x4 - x4 * 3.0 ) / std::sqrt(6);

    ASSERT_EQ ( result.size(), 1 );
    EXPECT_FLOAT_EQ ( result[0].real(), expected.real() );
    EXPECT_FLOAT_EQ ( result[0].imag(), expected.imag() );

}

TEST ( HermitePolynomials, H_2_1_1_2_0_3_2_ComplexInput ) {

    std::vector<size_t> indices { 2, 1, 1, 2, 0, 3 };

    std::complex<double> x1, x2;

    x1.real(  2.0 );
    x1.imag(  3.0 );

    x2.real( -5.0 );
    x2.imag(  2.0 );

    std::vector<std::complex<double>> X = { x1, x2 };

    size_t dim = 2;

    auto result =
        BasisFunctions::HermitePolynomials ( indices, X, dim );

    std::vector<std::complex<double>> expected;

    expected.push_back (
        x2 * ( x1 * x1 - 1.0 ) / std::sqrt(2)
    );

    expected.push_back (
        x1 * ( x2 * x2 - 1.0 ) / std::sqrt(2)
    );

    expected.push_back (
        ( x2 * x2 * x2 - x2 * 3.0 ) / std::sqrt(6)
    );

    ASSERT_EQ ( result.size(), expected.size() );

    for ( auto i = 0; i < expected.size(); i++ ) {
        
        EXPECT_FLOAT_EQ ( result[i].real(), expected[i].real() );
        EXPECT_FLOAT_EQ ( result[i].imag(), expected[i].imag() );

    }

}

TEST ( HermitePolynomials, H_Complex_ComplexInput ) {

    std::vector<size_t> indices { 2, 1, 1, 2, 0, 3 };

    std::complex<double> x1, x2, x3, x4;

    x1.real(  2.0 );
    x1.imag(  3.0 );

    x2.real( -5.0 );
    x2.imag(  2.0 );

    x3.real(  1.0 );
    x3.imag(  1.0 );

    x4.real( -1.0 );
    x4.imag(  2.0 );

    std::vector<std::complex<double>> X = { x1, x2, x3, x4 };

    size_t dim = 2;

    auto result =
        BasisFunctions::HermitePolynomials ( indices, X, dim );

    std::vector<std::complex<double>> expected;

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

    std::vector<size_t> indices { 2, 1, 1, 2, 0, 3 };

    std::complex<double> x1, x2;

    x1.real(  2.0 );
    x1.imag(  3.0 );

    x2.real( -5.0 );
    x2.imag(  2.0 );

    std::vector<std::complex<double>> X = { x1, x2 };

    size_t dim = 3;

    EXPECT_THROW ({

        try {

            auto result =
                BasisFunctions::HermitePolynomials ( indices, X, dim );

        } catch ( const std::exception& e ) {

            EXPECT_STREQ (
                "HermitePolynomials: num of samples not multiple of dimension",
                e.what()
            );

            throw;

        }


    }, std::runtime_error );


}

TEST ( HermitePolynomials, WrongIndicesSize ) {

    std::vector<size_t> indices { 2, 1, 1, 2, 0 };

    std::complex<double> x1, x2;

    x1.real(  2.0 );
    x1.imag(  3.0 );

    x2.real( -5.0 );
    x2.imag(  2.0 );

    std::vector<std::complex<double>> X = { x1, x2 };

    size_t dim = 2;

    EXPECT_THROW ({

        try {

            auto result =
                BasisFunctions::HermitePolynomials ( indices, X, dim );

        } catch ( const std::exception& e ) {

            EXPECT_STREQ (
                "HermitePolynomials: num of indices not multiple of dimension",
                e.what()
            );

            throw;

        }


    }, std::runtime_error );


}

TEST ( HermitePolynomials, ZeroDimension ) {

    std::vector<size_t> indices { 2, 1, 1, 2, 0, 3 };

    std::complex<double> x1, x2;

    x1.real(  2.0 );
    x1.imag(  3.0 );

    x2.real( -5.0 );
    x2.imag(  2.0 );

    std::vector<std::complex<double>> X = { x1, x2 };

    size_t dim = 0;

    EXPECT_THROW ({

        try {

            auto result =
                BasisFunctions::HermitePolynomials ( indices, X, dim );

        } catch ( const std::exception& e ) {

            EXPECT_STREQ (
                "HermitePolynomials: dimension must be positive",
                e.what()
            );

            throw;

        }


    }, std::runtime_error );


}

