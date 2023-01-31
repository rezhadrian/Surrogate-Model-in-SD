/**
  * @file testrunner.cpp
  *
  * @brief test of all non intrusive functionality 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include <gtest/gtest.h> 
#include "BasisFunctions.hpp" 
#include "MonteCarlo.hpp" 
#include "LinearAlgebra.hpp" 


TEST ( LA_MC, GenerateCorrelatedRVs ) {

    typedef linalg::DenseSymMatrix<double> SM;
    typedef linalg::DenseLTriangularMatrix<double> LT;
    typedef std::vector<double> Vector;

    size_t dim = 4;

    Vector StdNormRVs {

         0.1,  0.2,  0.3,  0.4, 
         0.4,  0.3,  0.2,  0.1, 
         0.5,  0.5,  0.4,  0.4, 
         0.3,  0.3,  0.2,  0.2, 
        -0.1, -0.2, -0.3, -0.4, 
        -0.4, -0.3, -0.2, -0.1, 
        -0.5, -0.5, -0.4, -0.4, 
        -0.3, -0.3, -0.2, -0.2 

    };


    SM Corr ( dim );

    Corr ( 0, 0 ) =  1.0000; Corr ( 1, 0 ) = -0.7601; 
    Corr ( 2, 0 ) =  0.9343; Corr ( 3, 0 ) =  0.9691; 
    Corr ( 1, 1 ) =  1.0000; Corr ( 2, 1 ) = -0.6293; 
    Corr ( 3, 1 ) = -0.7881; Corr ( 2, 2 ) =  1.0000; 
    Corr ( 3, 2 ) =  0.8920; Corr ( 3, 3 ) =  1.0000; 

    std::vector<std::function<double(double)>> ICDFs;


    ICDFs.emplace_back (

        [](const auto m) {
            return MonteCarlo::InvLogNormCDF<double> ( m, 0.0, 1.0 );
        }

    );

    ICDFs.emplace_back (

        [](const auto m) {
            return MonteCarlo::InvLogNormCDF<double> ( m, 0.1, 0.5 );
        }

    );

    ICDFs.emplace_back (

        [](const auto m) {
            return MonteCarlo::InvLogNormCDF<double> ( m, 0.3, 0.2 );
        }

    );

    ICDFs.emplace_back (

        [](const auto m) {
            return MonteCarlo::InvLogNormCDF<double> ( m, 0.3, 0.5 );
        }

    );

    auto L = linalg::Cholesky <LT,SM> ( Corr );

    
    auto result = 
        MonteCarlo::GenerateRVs<double,LT> ( StdNormRVs, L, ICDFs, dim );

    #ifdef MC_COMPLEX 
        std::vector<std::complex<double>> expected {
    #else
        Vector expected {
    #endif 

        1.105170918075648, 1.135389327255325, 1.410165929127596, 
        1.470494476447416, 1.491824697641270, 1.046499638624557, 
        1.485234831925537, 1.636457192059457, 1.648721270700128, 
        1.075113813834130, 1.541258060148288, 1.762968747246328, 
        1.349858807576003, 1.087037295008275, 1.457739423492359, 
        1.577350001291581, 0.904837418035960, 1.075756772448067, 
        1.292130778906117, 1.239119785606122, 0.670320046035639, 
        1.167131562286532, 1.226822022500117, 1.113453385295951, 
        0.606530659712634, 1.136068332899878, 1.182228237765192, 
        1.033551390650895, 0.740818220681718, 1.123607040686559, 
        1.249961941775021, 1.155177226930297
        
    };

    double tol = 0.000000000000001;

    ASSERT_EQ ( result.size(), expected.size() );

    for ( auto i = 0; i < result.size(); i++ ) {

        #ifdef MC_COMPLEX 
            EXPECT_NEAR ( result[i].real(), expected[i].real(), tol );
            EXPECT_NEAR ( result[i].imag(), expected[i].imag(), tol );
        #else 
            EXPECT_NEAR ( result[i], expected[i], tol );
        #endif

    }


}


// TEST ( Overall, FirstTry ) {
//
//     typedef std::complex<double> Complex; 
//     typedef std::vector<Complex> Vector;
//     typedef linalg::DenseMatrix<Complex> Matrix; 
//     typedef MonteCarlo::DenseSymMatrix <size_t, double> SM;
//     typedef MonteCarlo::DenseLTriangularMatrix <size_t,double> LT;
//
//
//
//     size_t dim  = 2;
//     size_t maxP = 1;
//
//     auto indices = BasisFunctions::MultiIndex ( dim, maxP );
//
//     size_t nSample = 9; 
//
//     auto LHSResult = 
//         MonteCarlo::LHS<
//
//             size_t, double, 
//             std::random_device, 
//             std::default_random_engine, 
//             std::uniform_real_distribution<double>, 
//             std::mt19937 
//
//         > ( nSample, dim );
//
//     ASSERT_EQ ( LHSResult.size() / dim, nSample );
//
//     MonteCarlo::ConvertLHStoStdNorm<size_t,double> ( LHSResult );
//
//     ASSERT_EQ ( LHSResult.size() / dim, nSample );
//
//     std::vector<std::function<double(double)>> ICDFs;
//
//     SM Corr ( dim );
//
//     Corr ( 0, 0 ) =  1.0000; Corr ( 1, 0 ) = -0.7601; 
//     Corr ( 1, 1 ) =  1.0000; 
//
//     ICDFs.emplace_back (
//
//         [](const auto m) {
//             return MonteCarlo::InvLogNormCDF<double> ( m, 0.0, 1.0 );
//         }
//
//     );
//
//     ICDFs.emplace_back (
//
//         [](const auto m) {
//             return MonteCarlo::InvLogNormCDF<double> ( m, 0.1, 0.5 );
//         }
//
//     );
//
//     std::vector<double> X = 
//         MonteCarlo::GenerateRVs < size_t, double,
//
//             SM, 
//             LT, 
//             std::function<double(double)>
//
//         > ( LHSResult, Corr, ICDFs, dim );
//
//     Vector XC ( X.size() );
//
//     std::transform ( 
//
//         X.begin(), X.end(), 
//         XC.begin(), 
//         [](const auto m ){
//             return Complex ( m );
//         }
//
//     );
//
//
//     Vector b {
//
//          1.0,  2.0,  3.0, 
//         -1.0, -3.0, -2.0, 
//          0.5,  0.1, -0.8 
//
//     };
//
//     auto data = BasisFunctions::HermitePolynomials ( indices, XC, dim );
//
//     auto nRow = X.size() / dim;
//     auto nCol = indices.size() / dim;
//
//     Matrix A ( nRow, nCol, data );
//
//     auto ATA = A.ConjTransProd ( A );
//     auto ATb = A.ConjTransProd ( b );
//
//     Vector x ( nCol, 0.0 );
//
//
//     auto Norm = [](const auto& v){return linalg::Euclidean<Vector,double>(v);};
//     using GS = linalg::GaussSeidel<Matrix,Vector,double,decltype(Norm)>;
//
//     double tol = 0.000000000001;
//     size_t MaxIter = 50;
//     GS GSSolver ( MaxIter, tol, Norm );
//
//     GSSolver.Solve ( ATA, ATb, x );
// }
//
