/**
  * @file LinAlg_test.cpp
  *
  * @brief test linear algebra functionalities 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "MonteCarlo.hpp" 
#include <gtest/gtest.h> 

// TEST ( DenseSymMatrix, LowerTriColumnAssignment ) {
//
//     typedef MonteCarlo::DenseSymMatrix<size_t, double> Matrix;
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
//     ASSERT_EQ ( A.data().size(), data.size() );
//
//     for ( auto i = 0; i < data.size(); i++ ) {
//
//         EXPECT_EQ ( A( Index1[i], Index2[i] ), data[i] );
//         EXPECT_EQ ( A( Index2[i], Index1[i] ), data[i] );
//
//     }
//
// }
//
// TEST ( DenseSymMatrix, UpperTriRowAssignment ) {
//
//     typedef MonteCarlo::DenseSymMatrix<size_t, double> Matrix;
//     typedef std::vector<double> Vector;
//
//     Matrix A ( 3 );
//
//     Vector Index1 { 0, 0, 0, 1, 1, 2 };
//     Vector Index2 { 0, 1, 2, 1, 2, 2 };
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
//     ASSERT_EQ ( A.data().size(), data.size() );
//
//     for ( auto i = 0; i < data.size(); i++ ) {
//
//         EXPECT_EQ ( A( Index1[i], Index2[i] ), data[i] );
//         EXPECT_EQ ( A( Index2[i], Index1[i] ), data[i] );
//
//     }
//
// }
//
// TEST ( DenseLTriangularMatrix, LowerTriColumnAssignment ) {
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
//
// TEST ( CholeskyDecomposition, DenseIdentity ) {
//
//     typedef MonteCarlo::DenseSymMatrix<size_t, double> Sym;
//     typedef MonteCarlo::DenseLTriangularMatrix<size_t, double> Tri;
//
//     size_t dim = 30;
//
//     Sym A ( dim );
//
//     for ( auto i = 0; i < dim; i++ ) {
//         A(i,i) = 1.0;
//     }
//
//     Tri L = MonteCarlo::CholeskyDecompose<Tri,Sym> ( A );
//
//     const Tri& l = L;
//
//     for ( auto i = 0; i < dim; i++ ) {
//     for ( auto j = 0; j < dim; j++ ) {
//
//         EXPECT_EQ ( A(i,j), l(i,j) );
//
//     }
//     }
//
// }
//
// TEST ( CholeskyDecomposition, DenseSPD ) {
//
//     typedef MonteCarlo::DenseSymMatrix<size_t, double> Sym;
//     typedef MonteCarlo::DenseLTriangularMatrix<size_t, double> Tri;
//
//     typedef std::vector<double> Vector;
//
//     size_t dim = 5;
//
//     Vector InputData {
//
//          9,  3,  0, 15,  0,
//          3,  5,  0,  7,  0,
//          0,  0, 16,  8,  0,
//         15,  7,  8, 39,  0,
//          0,  0,  0,  0,  1
//
//     };
//
//     ASSERT_EQ ( dim * dim, InputData.size() );
//
//     Sym A ( dim );
//
//     for ( auto i = 0; i < dim; i++ ) {
//     for ( auto j = 0; j < dim; j++ ) {
//
//         A(i,j) = InputData [i*dim + j];
//
//     }
//     }
//
//     Tri L = MonteCarlo::CholeskyDecompose<Tri,Sym> ( A );
//
//     Vector expected {
//
//         3, 0, 0, 0, 0,
//         1, 2, 0, 0, 0,
//         0, 0, 4, 0, 0,
//         5, 1, 2, 3, 0,
//         0, 0, 0, 0, 1
//
//     };
//
//     ASSERT_EQ ( dim * dim, expected.size() );
//
//     const Tri& l = L;
//
//     for ( auto i = 0; i < dim; i++ ) {
//     for ( auto j = 0; j < dim; j++ ) {
//
//         EXPECT_EQ ( l(i,j), expected[i*dim + j] );
//
//     }
//     }
//
// }
//
