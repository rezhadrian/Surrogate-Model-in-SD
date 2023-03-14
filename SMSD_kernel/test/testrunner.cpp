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


// TEST ( LA_MC, GenerateCorrelatedRVs ) {
//
//     typedef linalg::DenseSymMatrix<double> SM;
//     typedef linalg::DenseLTriangularMatrix<double> LT;
//     typedef std::vector<double> Vector;
//
//     size_t dim = 4;
//
//     Vector StdNormRVs {
//
//          0.1,  0.2,  0.3,  0.4, 
//          0.4,  0.3,  0.2,  0.1, 
//          0.5,  0.5,  0.4,  0.4, 
//          0.3,  0.3,  0.2,  0.2, 
//         -0.1, -0.2, -0.3, -0.4, 
//         -0.4, -0.3, -0.2, -0.1, 
//         -0.5, -0.5, -0.4, -0.4, 
//         -0.3, -0.3, -0.2, -0.2 
//
//     };
//
//
//     SM Corr ( dim );
//
//     Corr ( 0, 0 ) =  1.0000; Corr ( 1, 0 ) = -0.7601; 
//     Corr ( 2, 0 ) =  0.9343; Corr ( 3, 0 ) =  0.9691; 
//     Corr ( 1, 1 ) =  1.0000; Corr ( 2, 1 ) = -0.6293; 
//     Corr ( 3, 1 ) = -0.7881; Corr ( 2, 2 ) =  1.0000; 
//     Corr ( 3, 2 ) =  0.8920; Corr ( 3, 3 ) =  1.0000; 
//
//     std::vector<std::function<double(double)>> ICDFs;
//
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
//     ICDFs.emplace_back (
//
//         [](const auto m) {
//             return MonteCarlo::InvLogNormCDF<double> ( m, 0.3, 0.2 );
//         }
//
//     );
//
//     ICDFs.emplace_back (
//
//         [](const auto m) {
//             return MonteCarlo::InvLogNormCDF<double> ( m, 0.3, 0.5 );
//         }
//
//     );
//
//     auto L = linalg::Cholesky <LT,SM> ( Corr );
//
//     
//     auto result = 
//         MonteCarlo::GenerateRVs<double,LT> ( StdNormRVs, L, ICDFs, dim );
//
//     #ifdef MC_COMPLEX 
//         std::vector<std::complex<double>> expected {
//     #else
//         Vector expected {
//     #endif 
//
//         1.105170918075648, 1.135389327255325, 1.410165929127596, 
//         1.470494476447416, 1.491824697641270, 1.046499638624557, 
//         1.485234831925537, 1.636457192059457, 1.648721270700128, 
//         1.075113813834130, 1.541258060148288, 1.762968747246328, 
//         1.349858807576003, 1.087037295008275, 1.457739423492359, 
//         1.577350001291581, 0.904837418035960, 1.075756772448067, 
//         1.292130778906117, 1.239119785606122, 0.670320046035639, 
//         1.167131562286532, 1.226822022500117, 1.113453385295951, 
//         0.606530659712634, 1.136068332899878, 1.182228237765192, 
//         1.033551390650895, 0.740818220681718, 1.123607040686559, 
//         1.249961941775021, 1.155177226930297
//         
//     };
//
//     double tol = 0.000000000000001;
//
//     ASSERT_EQ ( result.size(), expected.size() );
//
//     for ( auto i = 0; i < result.size(); i++ ) {
//
//         #ifdef MC_COMPLEX 
//             EXPECT_NEAR ( result[i].real(), expected[i].real(), tol );
//             EXPECT_NEAR ( result[i].imag(), expected[i].imag(), tol );
//         #else 
//             EXPECT_NEAR ( result[i], expected[i], tol );
//         #endif
//
//     }
//
//
// }
//
// TEST ( Overall, FirstTry ) {
//
//     typedef std::complex<double> Complex; 
//     typedef std::vector<Complex> Vector;
//     typedef linalg::DenseMatrix<Complex> Matrix; 
//     typedef linalg::DenseSymMatrix <double> SM;
//     typedef linalg::DenseLTriangularMatrix <double> LT;
//
//     size_t dim = 2;
//     size_t nPoints = 25;
//
//     /** 
//       * Generate random variables: 
//       * 1. Perform sampling from uniform dist. (LHS) 
//       * 2. Convert LHS result to standard normal RVs
//       * 3. Decompose correlation matrix 
//       * 4. Use inverse transform to generate RVs 
//       */
//
//
//     // Perform sampling from uniform dist. (LHS) 
//     auto LHSResult = MonteCarlo::LHS<size_t,double> ( nPoints, dim );
//
//
//     // Convert LHS result to standard normal RVs 
//     MonteCarlo::ConvertLHStoStdNorm<size_t,double> ( LHSResult );
//
//
//     // Create correlation matrix and perform Cholesky decomposition 
//     SM Correlation ( dim );
//
//     Correlation ( 0, 0 ) =  1.0000; 
//     Correlation ( 1, 0 ) =  0.0000; Correlation ( 1, 1 ) = 1.0000;
//
//     LT L = linalg::Cholesky<LT> ( Correlation );
//
//
//     // Create a vector of ICDF functions for inverse transform 
//     std::vector< std::function < double(double) > > ICDFs;
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
//     auto RandomVariables = 
//         MonteCarlo::GenerateRVs<double,LT> ( LHSResult, L, ICDFs, dim );
//
//
//     /**
//       * Generate surrogate model
//       * 1. Evaluate basis functions using generated random variables 
//       * 2. Evaluate FRF model using generated random variables 
//       * 3. Perform least square fitting 
//       */
//
//
//     // Evaluate basis functions and assemble into matrix 
//     size_t maxP = 5;
//
//     auto indices = BasisFunctions::MultiIndex ( dim, maxP );
//
//     auto data = 
//         BasisFunctions::HermitePolynomials ( indices, RandomVariables, dim );
//
//     auto nRow = RandomVariables.size() / dim;
//     auto nCol = indices.size() / dim;
//
//     Matrix A ( nRow, nCol, data );
//
//
//     // Evaluate FRF model and assemble into vector 
//     auto model = []( const auto w, const auto& pars ) {
//         return MonteCarlo::EvaluateFRF<double, Complex> ( w, pars );
//     };
//
//     double omega = 10;
//
//     auto b = MonteCarlo::EvaluateModel<size_t, double, Complex, 
//          decltype(model)> ( omega, RandomVariables, model, dim );
//
//
//     // Perform least square fitting 
//     size_t MaxIter = 100;
//     double tol     = 1e-8;
//
//     auto ATA = A.ConjTransProd ( A ); 
//     auto ATb = A.ConjTransProd ( b );
//
//     Vector x ( nCol, 0.0 );
//
//     auto Norm = [](const auto& v){return linalg::Euclidean<Vector,double>(v);};
//
//     using GS = linalg::GaussSeidel<Matrix,Vector,double,decltype(Norm)>;
//     GS GSSolver ( MaxIter, tol, Norm );
//
//     GSSolver.Solve ( ATA, ATb, x );
//
//
// }
//
