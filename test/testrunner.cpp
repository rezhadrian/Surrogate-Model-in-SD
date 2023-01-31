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
