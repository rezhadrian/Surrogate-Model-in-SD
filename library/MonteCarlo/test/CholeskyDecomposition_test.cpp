/**
  * @file CholeskyDecomposition_test.cpp
  *
  * @brief test special matrix class for Monte Carlo  
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "MonteCarlo.hpp" 
#include <gtest/gtest.h>

TEST ( CholeskyDecomposition, DenseIdentity ) {

    typedef MonteCarlo::DenseSymMatrix<size_t, double> Sym;
    typedef MonteCarlo::DenseLTriangularMatrix<size_t, double> Tri;

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
    typedef MonteCarlo::DenseLTriangularMatrix<size_t, double> Tri;

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

