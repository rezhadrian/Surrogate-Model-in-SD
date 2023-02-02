/**
  * @file main.cpp
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "BasisFunctions.hpp" 
#include "MonteCarlo.hpp" 
// #include "LinearAlgebra.hpp" 
#include <Eigen/Dense> 

#include <iostream> 

using Eigen::MatrixXd;
using Eigen::VectorXcd;

int main () {
    //
    // typedef double Float;
    //
    // // typedef std::complex<Float> Complex; 
    // // typedef std::vector<Complex> Vector;
    // // typedef linalg::DenseMatrix<Complex> Matrix; 
    // // typedef linalg::DenseSymMatrix <Float> SM;
    // // typedef linalg::DenseLTriangularMatrix <Float> LT;
    //
    // size_t dim = 2;
    // size_t nPoints = 250000;
    //
    // /** 
    //   * Generate random variables: 
    //   * 1. Perform sampling from uniform dist. (LHS) 
    //   * 2. Convert LHS result to standard normal RVs
    //   * 3. Decompose correlation matrix 
    //   * 4. Use inverse transform to generate RVs 
    //   */
    //
    //
    // // Perform sampling from uniform dist. (LHS) 
    // auto LHSResult = MonteCarlo::LHS<size_t,Float> ( nPoints, dim );
    //
    //
    // // Convert LHS result to standard normal RVs 
    // MonteCarlo::ConvertLHStoStdNorm<size_t,Float> ( LHSResult );
    //
    //
    // // size_t maxP = 5;
    // // auto indices = BasisFunctions::MultiIndex ( dim, maxP );
    // //
    // // auto data = 
    // //     BasisFunctions::HermitePolynomials<size_t, Float, Float> (
    // //         indices, LHSResult, dim 
    // //     );
    // //
    // // std::cout << LHSResult.size() << std::endl << std::endl;
    //
    // // // Create correlation matrix and perform Cholesky decomposition 
    // // SM Correlation ( dim );
    // //
    // // Correlation ( 0, 0 ) =  1.0000; 
    // // Correlation ( 1, 0 ) =  0.0000; Correlation ( 1, 1 ) = 1.0000;
    // //
    // // LT L = linalg::Cholesky<LT> ( Correlation );
    // //
    // // auto NormRVs = 
    // //     MonteCarlo::CombineRVs<Float, LT> ( L, LHSResult, dim);
    // //
    // //
    // // // Create a vector of ICDF functions for inverse transform 
    // // std::vector< std::function < Float(Float) > > ICDFs;
    // //
    // // ICDFs.emplace_back ( 
    // //
    // //     [](const auto m) { 
    // //         return MonteCarlo::InvLogNormCDF<Float> ( 
    // //             m, 24.11948805, 0.099751 );
    // //     }
    // //
    // // );
    // //
    // // ICDFs.emplace_back ( 
    // //
    // //     [](const auto m) { 
    // //         return MonteCarlo::InvLogNormCDF<Float> ( 
    // //             m, 7.822797571, 0.0499687 
    // //         );
    // //     }
    // //
    // // );
    // //
    // // auto RandomVariables = 
    // //     MonteCarlo::GenerateRVs<Float,LT> ( LHSResult, L, ICDFs, dim );
    //
    //
    // /**
    //   * Generate surrogate model
    //   * 1. Evaluate basis functions using generated random variables 
    //   * 2. Evaluate FRF model using generated random variables 
    //   * 3. Perform least square fitting 
    //   */
    //
    //
    // // Evaluate basis functions and assemble into matrix 
    // //
    // // auto nRow = RandomVariables.size() / dim;
    // // auto nCol = indices.size() / dim;
    // //
    // //
    // //
    // //
    // //
    // //
    // // // Evaluate FRF model and assemble into vector 
    // // auto model = []( const auto w, const auto& pars ) {
    // //     return MonteCarlo::EvaluateFRF<Float, Complex> ( w, pars );
    // // };
    // //
    // // double omega = 32.044;
    // //
    // // auto b = MonteCarlo::EvaluateModel<size_t, Float, Complex, 
    // //      decltype(model)> ( omega, RandomVariables, model, dim );
    // //
    // //
    // // // Perform least square fitting 
    // // size_t MaxIter = 10000;
    // // double tol     = 1e-12;
    //
    // // auto ATA = A.ConjTransProd ( A ); 
    // // auto ATb = A.ConjTransProd ( b );
    // //
    // // Vector x ( nCol, 0.0 );
    // //
    // // auto Norm = 
    // //     [](const auto& v){return linalg::Euclidean<Vector,Float>(v);};
    // //
    // // using GS = linalg::GaussSeidel<Matrix,Vector,Float,decltype(Norm)>;
    // // GS GSSolver ( MaxIter, tol, Norm );
    // //
    // // GSSolver.Solve ( ATA, ATb, x );
    // //
    // // auto ToCompare = A * x;
    // //
    // //
    // // // for ( auto i = 0; i < 6; i++ ) {
    // // //
    // // //     std::cout << RandomVariables[i] << std::endl;
    // // //
    // // // }
    // //
    // // std::cout << std::endl << std::endl;
    // // std::cout << "Analytical model: ";
    // // std::cout << std::endl << std::endl;
    // //
    // // for ( auto i = 0; i < 5; i++ ) {
    // //
    // //     std::cout << b[i] << std::endl;
    // //
    // // }
    // //
    // // std::cout << std::endl << std::endl;
    // // std::cout << "Coefficients: ";
    // // std::cout << std::endl << std::endl;
    // //
    // // // output variables to console 
    // // for ( auto i = 0; i < x.size(); i++ ) {
    // //
    // //     std::cout << x[i] << std::endl;
    // //
    // // }
    // //
    // // std::cout << std::endl << std::endl;
    // // std::cout << "Surrogate model: ";
    // // std::cout << std::endl << std::endl;
    // //
    // // for ( auto i = 0; i < 5; i++ ) {
    // //
    // //     std::cout << ToCompare[i] << std::endl;
    // //
    // // }
    // //
    // // std::cout << std::endl;
    //

    return 0;
} // main 

