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
  * @file CholeskyDecomposition_test.cpp
  *
  * @brief test special matrix class for Monte Carlo  
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "SpecialMatrix.hpp" 
#include "CholeskyDecomposition.hpp" 
#include <gtest/gtest.h>

TEST ( CholeskyDecomposition, DenseIdentity ) {

    typedef MonteCarlo::DenseSymMatrix<size_t, double> Sym;
    typedef MonteCarlo::DenseTriangularMatrix<size_t, double> Tri;

    size_t dim = 30;

    Sym A ( dim );

    for ( auto i = 0; i < dim; i++ ) {
        A(i,i) = 1.0;
    }

    Tri L = MonteCarlo::CholeskyDecompose<Tri,Sym> ( A );

    const Tri& l = L;

    for ( auto i = 0; i < dim; i++ ) {
    for ( auto j = 0; j < dim; j++ ) {

        EXPECT_EQ ( A(i,j), l(i,j) );

    }
    }

}

TEST ( CholeskyDecomposition, DenseSPD ) {

    typedef MonteCarlo::DenseSymMatrix<size_t, double> Sym;
    typedef MonteCarlo::DenseTriangularMatrix<size_t, double> Tri;

    typedef std::vector<double> Vector;

    size_t dim = 5;

    Vector InputData {

         9,  3,  0, 15,  0,
         3,  5,  0,  7,  0,
         0,  0, 16,  8,  0,
        15,  7,  8, 39,  0,
         0,  0,  0,  0,  1

    };

    ASSERT_EQ ( dim * dim, InputData.size() );

    Sym A ( dim );

    for ( auto i = 0; i < dim; i++ ) {
    for ( auto j = 0; j < dim; j++ ) {

        A(i,j) = InputData [i*dim + j];

    }
    }

    Tri L = MonteCarlo::CholeskyDecompose<Tri,Sym> ( A );

    Vector expected {

        3, 0, 0, 0, 0,
        1, 2, 0, 0, 0, 
        0, 0, 4, 0, 0,
        5, 1, 2, 3, 0,
        0, 0, 0, 0, 1

    };

    ASSERT_EQ ( dim * dim, expected.size() );

    const Tri& l = L;

    for ( auto i = 0; i < dim; i++ ) {
    for ( auto j = 0; j < dim; j++ ) {

        EXPECT_EQ ( l(i,j), expected[i*dim + j] );

    }
    }

}

