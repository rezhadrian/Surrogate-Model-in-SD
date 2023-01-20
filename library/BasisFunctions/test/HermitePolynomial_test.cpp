/**
  * This file is part of a Surrogate Model library for structural dynamics. 
  *
  * The library is free: you can distribute it and/or modify it.
  * The library is distributed in the hope that it can be useful and helpful,
  * particularly for students and self - learning individuals.
  * 
  * The library comes WITHOUT ANY WARRANTY: without even the implied warranty 
  * of MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE.
  *
  * The library is licensed under the GNU General Public License v3.0
  *
  *
  * The library is based on works by:
  * - Felix Schneider
  * - Iason Papaioannou
  * 
  * Chair of Structural Mechanics 
  * Engineering Risk Analysis Group 
  *
  * Technische Universitaet Muenchen
  *
  * www.cee.ed.tum.de/bm
  * www.cee.ed.tum.de/era 
  */ 

/**
  * @file HermitePolynomial_test.cpp
  *
  * @brief test output of non-normalized HermitePolynomial function 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "HermitePolynomials.hpp" 
#include <gtest/gtest.h> 

TEST ( HermitePolynomial, H_0_ComplexInput ) {

    size_t index = 0;

    std::complex<double> x;
    x.real (  3.0 );
    x.imag (  2.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    std::complex<double> expected;
    expected.real ( 1.0 );
    expected.imag ( 0.0 );

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_0_RealInput ) {

    size_t index = 0;

    std::complex<double> x;
    x.real (  7.0 );
    x.imag (  0.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    std::complex<double> expected;
    expected.real ( 1.0 );
    expected.imag ( 0.0 );

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_0_ImagInput ) {

    size_t index = 0;

    std::complex<double> x;
    x.real (  0.0 );
    x.imag (  9.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    std::complex<double> expected;
    expected.real ( 1.0 );
    expected.imag ( 0.0 );

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_0_NegInput ) {

    size_t index = 0;

    std::complex<double> x;
    x.real ( -1.0 );
    x.imag (  0.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    std::complex<double> expected;
    expected.real ( 1.0 );
    expected.imag ( 0.0 );

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_0_ZeroInput ) {

    size_t index = 0;

    std::complex<double> x;
    x.real (  0.0 );
    x.imag (  0.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    std::complex<double> expected;
    expected.real ( 1.0 );
    expected.imag ( 0.0 );

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_1_ComplexInput ) {

    size_t index = 1;

    std::complex<double> x;
    x.real (  2.0 );
    x.imag (  2.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    auto expected = x;

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_1_RealInput ) {

    size_t index = 1;

    std::complex<double> x;
    x.real (  8.0 );
    x.imag (  0.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    auto expected = x;

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_1_ImagInput ) {

    size_t index = 1;

    std::complex<double> x;
    x.real (  0.0 );
    x.imag (  2.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    auto expected = x;

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_1_NegInput ) {

    size_t index = 1;

    std::complex<double> x;
    x.real (  0.0 );
    x.imag ( -5.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    auto expected = x;

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_1_ZeroInput ) {

    size_t index = 1;

    std::complex<double> x;
    x.real (  0.0 );
    x.imag (  0.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    auto expected = x;

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_2_ComplexInput ) {

    size_t index = 2;

    std::complex<double> x;
    x.real (  1.0 );
    x.imag (  3.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    auto expected = x * x - 1.0;

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_3_ComplexInput ) {

    size_t index = 3;

    std::complex<double> x;
    x.real (  4.0 );
    x.imag (  5.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    auto expected = x * x * x - x * 3.0;

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

TEST ( HermitePolynomial, H_4_ComplexInput ) {

    size_t index = 4;

    std::complex<double> x;
    x.real ( -9.0 );
    x.imag (  2.0 );

    auto result = BasisFunctions::HermitePolynomial ( index, x );

    auto expected = x * x * x * x - x * x * 6.0 + 3.0;

    EXPECT_FLOAT_EQ ( result.real(), expected.real() );
    EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

}

