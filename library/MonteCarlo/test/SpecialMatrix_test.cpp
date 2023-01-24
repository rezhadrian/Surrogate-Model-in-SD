/**
  * @file SpecialMatrix_test.cpp
  *
  * @brief test special matrix class for Monte Carlo  
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "MonteCarlo.hpp" 
#include <gtest/gtest.h> 

TEST ( DenseSymMatrix, LowerTriColumnAssignment ) {

    typedef MonteCarlo::DenseSymMatrix<size_t, double> Matrix;
    typedef std::vector<double> Vector;

    Matrix A ( 3 );

    Vector Index1 { 0, 1, 2, 1, 2, 2 };
    Vector Index2 { 0, 0, 0, 1, 1, 2 };

    Vector data { 1, 2, 3, 4, 5, 6 };


    for ( auto i = 0; i < data.size(); i++ ) {

        A( Index1[i], Index2[i] ) = data[i];

    }

    ASSERT_EQ ( A.data().size(), data.size() );

    for ( auto i = 0; i < data.size(); i++ ) {

        EXPECT_EQ ( A( Index1[i], Index2[i] ), data[i] );
        EXPECT_EQ ( A( Index2[i], Index1[i] ), data[i] );

    }

}

TEST ( DenseSymMatrix, UpperTriRowAssignment ) {

    typedef MonteCarlo::DenseSymMatrix<size_t, double> Matrix;
    typedef std::vector<double> Vector;

    Matrix A ( 3 );

    Vector Index1 { 0, 0, 0, 1, 1, 2 };
    Vector Index2 { 0, 1, 2, 1, 2, 2 };

    Vector data { 1, 2, 3, 4, 5, 6 };


    for ( auto i = 0; i < data.size(); i++ ) {

        A( Index1[i], Index2[i] ) = data[i];

    }

    ASSERT_EQ ( A.data().size(), data.size() );

    for ( auto i = 0; i < data.size(); i++ ) {

        EXPECT_EQ ( A( Index1[i], Index2[i] ), data[i] );
        EXPECT_EQ ( A( Index2[i], Index1[i] ), data[i] );

    }

}

TEST ( DenseUTriangularMatrix, UpperTriRowAssignment ) {

    typedef MonteCarlo::DenseUTriangularMatrix<size_t, double> Matrix;
    typedef std::vector<double> Vector;

    Matrix A ( 3 );

    Vector Index1 { 0, 0, 0, 1, 1, 2 };
    Vector Index2 { 0, 1, 2, 1, 2, 2 };

    Vector data { 1, 2, 3, 4, 5, 6 };


    for ( auto i = 0; i < data.size(); i++ ) {

        A( Index1[i], Index2[i] ) = data[i];

    }

    const Matrix& a = A;

    Vector LowerRow {  1, 2, 2 };
    Vector LowerCol {  0, 0, 1 };

    for ( auto i = 0; i < data.size(); i++ ) {

        EXPECT_EQ ( A( Index1[i], Index2[i] ), data[i] );

    }

    for ( auto i = 0; i < LowerRow.size(); i++ ) {

        EXPECT_EQ ( a( LowerRow[i], LowerCol[i] ), 0.0 );

    }

}

// TEST ( DenseTriangularMatrix, LowerTriColumnAssignment ) {
//
//     typedef MonteCarlo::DenseLTriangularMatrix<size_t, double> Matrix;
//     typedef std::vector<double> Vector;
//
//     Matrix A ( 3 );
//
//     Vector Index1 { 0, 1, 2, 1, 2, 2 };
//     Vector Index2 { 0, 0, 0, 1, 1, 2 };
//
//     Vector data { 1, 2, 3, 4, 5, 6 };
//
//
//     for ( auto i = 0; i < data.size(); i++ ) {
//
//         A( Index1[i], Index2[i] ) = data[i];
//
//     }
//
//     const Matrix& a = A;
//
//     Vector UpperRow {  0, 0, 1 };
//     Vector UpperCol {  1, 2, 2 };
//
//     for ( auto i = 0; i < data.size(); i++ ) {
//
//         EXPECT_EQ ( A( Index1[i], Index2[i] ), data[i] );
//
//     }
//
//     for ( auto i = 0; i < UpperCol.size(); i++ ) {
//
//         EXPECT_EQ ( a( UpperRow[i], UpperCol[i] ), 0.0 );
//
//     }
//
// }

