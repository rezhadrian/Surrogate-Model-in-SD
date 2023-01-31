/**
  * @file GaussSeidel_test.cpp
  *
  * @brief test Gauss Seidel iterative solver 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "LinearAlgebra.hpp" 
#include <gtest/gtest.h> 

TEST ( GaussSeidel, RealSysOfEqs ) {

    #ifdef LA_COMPLEX 
        typedef std::vector<std::complex<double>> Vector;
        typedef linalg::DenseMatrix<std::complex<double>> Matrix; 
    #else
        typedef std::vector<double> Vector;
        typedef linalg::DenseMatrix<double> Matrix; 
    #endif

    size_t nRow = 3;
    size_t nCol = 3;

    Vector data {

        15.0,  2.0,  3.0, 
         4.0, 11.0,  6.0, 
         1.0,  2.0,  8.0 

    };

    Matrix A ( nRow, nCol, data );

    Vector b { 1.0, 1.0, 1.0 };
    Vector x ( nCol, 0.0 );

    auto Norm = [](const auto& v){return linalg::Euclidean<Vector,double>(v);};

    using GS = linalg::GaussSeidel<Matrix,Vector,double,decltype(Norm)>;

    double tol = 0.00000001;
    size_t MaxIter = 50;

    GS GSSolver ( MaxIter, tol, Norm );

    GSSolver.Solve ( A, b, x );

    Vector expected { 

        0.041705282669138, 
        0.012048192771084, 
        0.116774791473587

    };

    ASSERT_EQ ( x.size(), expected.size() );

    for ( auto i = 0; i < x.size(); i++ ) {

        #ifdef LA_COMPLEX 
            EXPECT_NEAR ( x[i].real(), expected[i].real(), tol );
            EXPECT_NEAR ( x[i].imag(), expected[i].imag(), tol );
        #else 
            EXPECT_NEAR ( x[i], expected[i], tol );
        #endif

    }

}

#ifdef LA_COMPLEX 
TEST ( GaussSeidel, ComplexSysOfEqs ) {

    typedef std::complex<double> Complex;
    typedef std::vector<Complex> Vector;
    typedef linalg::DenseMatrix<Complex> Matrix; 

    size_t nRow = 3;
    size_t nCol = 3;

    Vector data {

        Complex ( 0.0, 7.0 ), Complex ( 1.0,-3.0), Complex ( 1.0, 0.0 ), 
        Complex (-1.0, 0.0 ), Complex ( 8.0, 3.0), Complex (-4.0, 0.0 ), 
        Complex ( 0.0, 0.0 ), Complex ( 0.0, 0.0), Complex ( 5.0, 0.0 ) 

    };

    Matrix A ( nRow, nCol, data );

    Vector b { 1.0, 1.0, 1.0 };
    Vector x ( nCol, 0.0 );

    auto Norm = [](const auto& v){return linalg::Euclidean<Vector,double>(v);};

    using GS = linalg::GaussSeidel<Matrix,Vector,double,decltype(Norm)>;

    double tol = 0.00000001;
    size_t MaxIter = 50;

    GS GSSolver ( MaxIter, tol, Norm );

    GSSolver.Solve ( A, b, x );

    Vector expected { 

        Complex ( 0.100155811779371, -0.124587098784668 ), 
        Complex ( 0.203116235587410, -0.091741975693362 ), 
        Complex ( 0.200000000000000,  0.000000000000000 )

    };

    ASSERT_EQ ( x.size(), expected.size() );

    for ( auto i = 0; i < x.size(); i++ ) {

        EXPECT_NEAR ( x[i].real(), expected[i].real(), tol );
        EXPECT_NEAR ( x[i].imag(), expected[i].imag(), tol );

    }

}
#endif

