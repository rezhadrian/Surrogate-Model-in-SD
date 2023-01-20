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

    auto result =
        BasisFunctions::HermitePolynomial<size_t, std::complex<double>> (
            0, std::complex<double> (1,1)
    );

    EXPECT_FLOAT_EQ ( result.real(), 1.0 );
    EXPECT_FLOAT_EQ ( result.imag(), 0.0 );

}

TEST ( HermitePolynomial, H_0_RealInput ) {

    auto result =
        BasisFunctions::HermitePolynomial<size_t, std::complex<double>> (
            0, std::complex<double> (1,0)
    );

    EXPECT_FLOAT_EQ ( result.real(), 1.0 );
    EXPECT_FLOAT_EQ ( result.imag(), 0.0 );

}

TEST ( HermitePolynomial, H_0_ImagInput ) {

    auto result =
        BasisFunctions::HermitePolynomial<size_t, std::complex<double>> (
            0, std::complex<double> (0,1)
    );

    EXPECT_FLOAT_EQ ( result.real(), 1.0 );
    EXPECT_FLOAT_EQ ( result.imag(), 0.0 );

}

TEST ( HermitePolynomial, H_0_NegInput ) {

    auto result =
        BasisFunctions::HermitePolynomial<size_t, std::complex<double>> (
            0, std::complex<double> (0,-1)
    );

    EXPECT_FLOAT_EQ ( result.real(), 1.0 );
    EXPECT_FLOAT_EQ ( result.imag(), 0.0 );

}

TEST ( HermitePolynomial, H_0_ZeroInput ) {

    auto result =
        BasisFunctions::HermitePolynomial<size_t, std::complex<double>> (
            0, std::complex<double> (0,0)
    );

    EXPECT_FLOAT_EQ ( result.real(), 1.0 );
    EXPECT_FLOAT_EQ ( result.imag(), 0.0 );

}

TEST ( HermitePolynomial, H_1_ComplexInput ) {

    auto result =
        BasisFunctions::HermitePolynomial<size_t, std::complex<double>> (
            1, std::complex<double> (1,1)
    );

    EXPECT_FLOAT_EQ ( result.real(), 1.0 );
    EXPECT_FLOAT_EQ ( result.imag(), 1.0 );

}

TEST ( HermitePolynomial, H_1_RealInput ) {

    auto result =
        BasisFunctions::HermitePolynomial<size_t, std::complex<double>> (
            1, std::complex<double> (3,0)
    );

    EXPECT_FLOAT_EQ ( result.real(), 3.0 );
    EXPECT_FLOAT_EQ ( result.imag(), 0.0 );

}

TEST ( HermitePolynomial, H_1_ImagInput ) {

    auto result =
        BasisFunctions::HermitePolynomial<size_t, std::complex<double>> (
            1, std::complex<double> (0,5)
    );

    EXPECT_FLOAT_EQ ( result.real(), 0.0 );
    EXPECT_FLOAT_EQ ( result.imag(), 5.0 );

}

TEST ( HermitePolynomial, H_1_NegInput ) {

    auto result =
        BasisFunctions::HermitePolynomial<size_t, std::complex<double>> (
            1, std::complex<double> (0,-5)
    );

    EXPECT_FLOAT_EQ ( result.real(),  0.0 );
    EXPECT_FLOAT_EQ ( result.imag(), -5.0 );

}

TEST ( HermitePolynomial, H_1_ZeroInput ) {

    auto result =
        BasisFunctions::HermitePolynomial<size_t, std::complex<double>> (
            1, std::complex<double> (0,0)
    );

    EXPECT_FLOAT_EQ ( result.real(), 0.0 );
    EXPECT_FLOAT_EQ ( result.imag(), 0.0 );

}

