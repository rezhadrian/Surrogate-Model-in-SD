/**
  * @file Cholesky_test.cpp
  *
  * @brief test implementations of Cholesky decomposition 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "LinearAlgebra.hpp" 
#include <gtest/gtest.h> 

TEST ( CholeskyDecomposition, DenseIdentity ) {

    typedef linalg::DenseSymMatrix<double> Sym;
    typedef linalg::DenseLTriangularMatrix<double> Tri;

    size_t dim = 30;

    Sym A ( dim );

    for ( auto i = 0; i < dim; i++ ) {
        A(i,i) = 1.0;
    }

    Tri L = linalg::Cholesky<Tri,Sym> ( A );

    const Tri& l = L;

    for ( auto i = 0; i < dim; i++ ) {
    for ( auto j = 0; j < dim; j++ ) {

        EXPECT_EQ ( A(i,j), l(i,j) );

    }
    }

}

TEST ( CholeskyDecomposition, DenseSPD ) {

    typedef linalg::DenseSymMatrix<double> Sym;
    typedef linalg::DenseLTriangularMatrix<double> Tri;

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

    Tri L = linalg::Cholesky<Tri,Sym> ( A );

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

TEST ( CholeskyDecomposition, NonSquareMatrix ) {

    typedef linalg::DenseMatrix<double> DM;
    typedef linalg::DenseLTriangularMatrix<double> LT;

    DM A ( 30, 40 );

    try {

        auto result = linalg::Cholesky<LT,DM> ( A );

    } catch ( const std::exception& e ) {

        EXPECT_STREQ (
            "Cholesky: Matrix is not square", 
            e.what()
        );

    }

}
TEST ( CholeskyDecomposition, NonSPDMatrix ) {

    typedef linalg::DenseMatrix<double> DM;
    typedef linalg::DenseLTriangularMatrix<double> LT;

    DM A ( 30, 30 );

    try {

        auto result = linalg::Cholesky<LT,DM> ( A );

    } catch ( const std::exception& e ) {

        EXPECT_STREQ (
            "Cholesky: input matrix is not pos. definite", 
            e.what()
        );

    }

}

