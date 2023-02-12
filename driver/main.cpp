/**
  * @file main.cpp
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "BasisFunctions.hpp" 
#include "MonteCarlo.hpp" 
#include "SurrogateModel.hpp" 
#include <Eigen/Dense> 

#include <iostream> 

using Eigen::MatrixXd;
using Eigen::VectorXcd;

int main () {

    typedef double Float; 
    typedef std::complex<Float> Complex; 

    typedef std::vector<Float>   VectorF;
    typedef std::vector<Complex> VectorC;

    typedef Eigen::Vector<Complex, Eigen::Dynamic> VectorXC;

    typedef Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> MatrixXF;
    typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> MatrixXC;

    typedef Eigen::DiagonalMatrix<Complex, Eigen::Dynamic> DiagMatrixXC;

    typedef Analytical::MassSpringDamper<size_t,Float,Complex> Model;
    typedef Surrogate::NonIntrusivePCE SuModel;

    size_t dim = 2;
    size_t nPoints = 25;

    VectorF Masses { 2.0, 1.0 };
    VectorF Damper { 2.0, 1.0 };
    VectorF Spring { 2.0, 1.0 };

    VectorC Force  { 1.0, 0.0 };

    VectorF AdditionalSpring { 0.3, 0.5, -0.1, 0.2 };
    VectorC AddSpring { 0.3, 0.5, -0.1, 0.2 };

    Float Omega = 5.0; 

    Model SDModel ( Masses, Damper, Spring );

    SuModel PCE ( &SDModel, Omega );

    PCE.SetIndices (
        5, 3
    );

    PCE.Train ( Force, AdditionalSpring );

    PCE.PrintCoeffs ();

    auto dispResult = PCE.ComputeResponse ( AddSpring );

    for ( auto i = 0; i < dispResult.size(); i++ ) {

        std::cout << dispResult[i] << std::endl;

    }

    // std::cout << 
    //     "=========================================\n\n";


    // /** 
    //   * Generate random variables: 
    //   * 1. Perform sampling from uniform dist. (LHS) 
    //   * 2. Convert LHS result to standard normal RVs
    //   * 3. Decompose correlation matrix 
    //   * 4. Use inverse transform to generate RVs 
    //   */
    //
    // // Perform sampling from uniform dist. (LHS) 
    // auto LHSResult = MonteCarlo::LHS<size_t,Float> ( nPoints, dim );
    //
    //
    // // Convert LHS result to standard normal RVs 
    // MonteCarlo::ConvertLHStoStdNorm<size_t,Float> ( LHSResult );
    //
    //
    //
    // MatrixXF Correl ( 2, 2 );
    //
    // Correl ( 0, 0 ) = 1.0;
    // Correl ( 1, 1 ) = 1.0;
    //
    //
    // auto X = MonteCarlo::CombineRVs<Float,Complex> ( 
    //
    //         Correl, 
    //         LHSResult, 
    //         dim 
    //
    // );
    //
    // // Create a vector of ICDF functions for inverse transform 
    // std::vector< std::function < Float(Float) > > ICDFs;
    //
    // boost::math::lognormal dist1 ( Float(24.11948805), Float(0.099751) );
    //
    // ICDFs.emplace_back ( 
    //
    //     [dist1](const auto m) { 
    //         return quantile ( dist1, m );
    //     }
    //
    // );
    //
    // boost::math::lognormal dist2 ( Float(7.822797571), Float(0.0499687) );
    //
    // ICDFs.emplace_back ( 
    //
    //     [dist2](const auto m) { 
    //         return quantile ( dist2, m );
    //     }
    //
    // );
    //
    // auto eX = MonteCarlo::GenerateRVs<size_t,Float> ( 
    //
    //     LHSResult, 
    //     Correl, 
    //     ICDFs, 
    //     dim 
    //
    // );
    //
    // /**
    //   * Generate surrogate model
    //   * 1. Evaluate basis functions using generated random variables 
    //   * 2. Evaluate FRF model using generated random variables 
    //   * 3. Perform least square fitting 
    //   */
    // 
    // size_t Mp = 2;
    // size_t Mq = 2;
    //
    // auto indicesP = BasisFunctions::MultiIndex ( dim, Mp );
    // auto indicesQ = BasisFunctions::MultiIndex ( dim, Mq );
    //
    // auto nP = indicesP.size() / dim;
    // auto nQ = indicesQ.size() / dim;
    //
    // // Evaluate basis functions 
    // auto psiP = 
    //     BasisFunctions::HermitePolynomials<size_t, Float, Complex> (
    //         indicesP, X, dim 
    //     );
    //
    // auto psiQ = 
    //     BasisFunctions::HermitePolynomials<size_t, Float, Complex> (
    //         indicesQ, X, dim 
    //     );
    //
    // // Evaluate FRF model and assemble into vector 
    // auto model = []( const auto w, const auto& pars ) {
    //     return MonteCarlo::EvaluateFRF<Float, Complex> ( w, pars );
    // };
    //
    // Float omega = 30;
    //
    // auto m = MonteCarlo::EvaluateModel <
    //
    //     size_t,Float,Complex, decltype(model)
    //
    // > (
    //
    //     omega, 
    //     eX, 
    //     model, 
    //     dim 
    //
    // );
    //
    // Eigen::Map<MatrixXC> PsiP ( psiP.data(), nP, nPoints );
    // Eigen::Map<MatrixXC> PsiQ ( psiQ.data(), nQ, nPoints );
    //
    // Eigen::Map<VectorXC> mv ( m.data(), nPoints );
    // DiagMatrixXC M = mv.asDiagonal();
    // DiagMatrixXC Mh = mv.conjugate().asDiagonal();
    //
    // MatrixXC Aupper ( nP, nP + nQ );
    // Aupper << ( PsiP * PsiP.transpose() ), ( -PsiP * M * PsiQ.transpose() );
    //
    // MatrixXC Alower ( nQ, nP + nQ );
    // Alower << ( -PsiQ * Mh * PsiP.transpose() ), ( PsiQ * Mh * M * PsiP.transpose() );
    //
    //
    // MatrixXC A ( nP + nQ, nP + nQ );
    // A << Aupper, Alower;
    //
    // MatrixXC V = A.bdcSvd( Eigen::ComputeFullV ).matrixV();
    //
    // // VectorXC r = V.col( V.cols() - 1 );
    //
    // std::cout << V.col( V.cols() - 1 ) << std::endl;

}

// int main () {
//
//     typedef double Float; 
//     typedef std::complex<Float> Complex; 
//
//     typedef std::vector<Float>   VectorF;
//     typedef std::vector<Complex> VectorC;
//
//     typedef Eigen::Vector<Complex, Eigen::Dynamic> VectorXC;
//
//     typedef Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> MatrixXF;
//     typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> MatrixXC;
//
//     size_t dim = 2; 
//     size_t nPoints = 25;
//
//
//     /** 
//       * Generate random variables: 
//       * 1. Perform sampling from uniform dist. (LHS) 
//       * 2. Convert LHS result to standard normal RVs
//       * 3. Decompose correlation matrix 
//       * 4. Use inverse transform to generate RVs 
//       */
//
//     // Perform sampling from uniform dist. (LHS) 
//
//     auto LHSResult = MonteCarlo::LHS<size_t,Float> ( nPoints, dim );
//
//
//     // Convert LHS result to standard normal RVs 
//     MonteCarlo::ConvertLHStoStdNorm<size_t,Float> ( LHSResult );
//
//
//     size_t maxP = 8;
//     auto indices = BasisFunctions::MultiIndex ( dim, maxP );
//
//     MatrixXF Correl ( 2, 2 );
//
//     Correl ( 0, 0 ) = 1.0;
//     Correl ( 1, 1 ) = 1.0;
//
//
//     auto StdNormRVs = MonteCarlo::CombineRVs<Float,Complex> ( 
//
//             Correl, 
//             LHSResult, 
//             dim 
//
//     );
//
//
//
//     // Create a vector of ICDF functions for inverse transform 
//     std::vector< std::function < Float(Float) > > ICDFs;
//
//     boost::math::lognormal dist1 ( Float(24.11948805), Float(0.099751) );
//
//     ICDFs.emplace_back ( 
//
//         [dist1](const auto m) { 
//             return quantile ( dist1, m );
//         }
//
//     );
//
//     boost::math::lognormal dist2 ( Float(7.822797571), Float(0.0499687) );
//
//     ICDFs.emplace_back ( 
//
//         [dist2](const auto m) { 
//             return quantile ( dist2, m );
//         }
//
//     );
//
//     auto RVs = MonteCarlo::GenerateRVs<size_t,Float> ( 
//
//         LHSResult, 
//         Correl, 
//         ICDFs, 
//         dim 
//
//     );
//
//
//     /**
//       * Generate surrogate model
//       * 1. Evaluate basis functions using generated random variables 
//       * 2. Evaluate FRF model using generated random variables 
//       * 3. Perform least square fitting 
//       */
//
//     // Evaluate basis functions 
//     auto Basis = 
//         BasisFunctions::HermitePolynomials<size_t, Float, Complex> (
//             indices, StdNormRVs, dim 
//         );
//
//
//     // Evaluate FRF model and assemble into vector 
//     auto model = []( const auto w, const auto& pars ) {
//         return MonteCarlo::EvaluateFRF<Float, Complex> ( w, pars );
//     };
//
//     Float omega = 30;
//
//     auto bStdVec = MonteCarlo::EvaluateModel <
//
//         size_t,Float,Complex, decltype(model)
//
//     > (
//
//         omega, 
//         RVs, 
//         model, 
//         dim 
//
//     );
//
//     auto nRow = bStdVec.size();
//     auto nCol = indices.size() / dim;
//
//     Eigen::Map<MatrixXC> A ( Basis.data(), nCol, nRow );
//
//     Eigen::Map<VectorXC> b ( bStdVec.data(), bStdVec.size() );
//
//     VectorXC x = ( A * A.transpose() ).ldlt().solve( A * b);
//
//     std::cout << x << std::endl;
//
//
//     return 0;
// } // main 
//
