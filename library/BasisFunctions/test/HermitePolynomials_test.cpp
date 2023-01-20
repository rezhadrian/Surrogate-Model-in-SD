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
  * @file HermitePolynomials_test.cpp
  *
  * @brief test output of HermitePolynomials function 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "HermitePolynomials.hpp" 
#include <gtest/gtest.h> 

TEST ( HermitePolynomials, H_0_0_2_ZeroInput ) {

    std::vector<size_t> indices {0,0};
    std::vector<std::complex<double>> X;

    X.push_back ( std::complex<double> (0,0) );
    X.push_back ( std::complex<double> (0,0) );

    size_t dim = 2;

    auto result =
        BasisFunctions::HermitePolynomials ( indices, X, dim );

    EXPECT_EQ ( result.size(), 1 );
    EXPECT_FLOAT_EQ ( result[0].real(), 1.0 );
    EXPECT_FLOAT_EQ ( result[0].imag(), 0.0 );

}

TEST ( HermitePolynomials, H_1_1_2_ZeroInput ) {

    std::vector<size_t> indices {1,1};
    std::vector<std::complex<double>> X;

    X.push_back ( std::complex<double> (0,0) );
    X.push_back ( std::complex<double> (0,0) );

    size_t dim = 2;

    auto result =
        BasisFunctions::HermitePolynomials ( indices, X, dim );

    EXPECT_EQ ( result.size(), 1 );
    EXPECT_FLOAT_EQ ( result[0].real(), 0.0 );
    EXPECT_FLOAT_EQ ( result[0].imag(), 0.0 );

}

TEST ( HermitePolynomials, H_0_1_2_ZeroInput ) {

    std::vector<size_t> indices {0,1};
    std::vector<std::complex<double>> X;

    X.push_back ( std::complex<double> (0,0) );
    X.push_back ( std::complex<double> (0,0) );

    size_t dim = 2;

    auto result =
        BasisFunctions::HermitePolynomials ( indices, X, dim );

    EXPECT_EQ ( result.size(), 1 );
    EXPECT_FLOAT_EQ ( result[0].real(), 0.0 );
    EXPECT_FLOAT_EQ ( result[0].imag(), 0.0 );

}

