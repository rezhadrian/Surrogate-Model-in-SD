/**
  * @file LeastSquare_test.cpp
  *
  * @brief test functions to solve least square problems 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "LinearAlgebra.hpp" 
#include <gtest/gtest.h> 

TEST ( LeastSquare, OverDetermined ) {

    typedef std::complex<double> Complex;
    typedef std::vector<Complex> Vector;
    typedef linalg::DenseMatrix<Complex> Matrix;

    Vector data {

        Complex ( 3.0, 0.0 ), Complex ( 0.0, 0.0 ), 
        Complex ( 0.0, 0.0 ), Complex ( 3.0, 0.0 ), 
        Complex ( 1.0, 0.0 ), Complex ( 1.0, 0.0 )

    };

    Matrix A ( 3, 2, data );

    Vector b {

        Complex ( 1.0, 0.0 ), 
        Complex ( 0.0, 1.0 ), 
        Complex ( 3.0, 0.0 )

    };

    auto ATA = A.ConjTransProd ( A );
    auto ATb = A.ConjTransProd ( b );

    Vector x ( 2, 0.0 );

    auto Norm = [](const auto& v){return linalg::Euclidean<Vector,double>(v);};
    using GS = linalg::GaussSeidel<Matrix,Vector,double,decltype(Norm)>;

    double tol = 0.000000000001;
    size_t MaxIter = 50;

    GS GSSolver ( MaxIter, tol, Norm );

    GSSolver.Solve ( ATA, ATb, x );

    Vector expected {

        Complex ( 0.575757575757576, -0.030303030303030 ),
        Complex ( 0.242424242424242,  0.303030303030303 )

    };

    ASSERT_EQ ( x.size(), expected.size() );

    for ( auto i = 0; i < x.size(); i++ ) {

        EXPECT_NEAR ( x[i].real(), expected[i].real(), tol );
        EXPECT_NEAR ( x[i].imag(), expected[i].imag(), tol );

    }

} 

