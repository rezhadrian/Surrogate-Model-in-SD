/**
  * @file main.cpp
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "BasisFunctions.hpp" 
#include "MonteCarlo.hpp" 
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

    size_t dim = 2; 
    size_t nPoints = 2500;


    /** 
      * Generate random variables: 
      * 1. Perform sampling from uniform dist. (LHS) 
      * 2. Convert LHS result to standard normal RVs
      * 3. Decompose correlation matrix 
      * 4. Use inverse transform to generate RVs 
      */

    // Perform sampling from uniform dist. (LHS) 

    auto LHSResult = MonteCarlo::LHS<size_t,Float> ( nPoints, dim );


    // Convert LHS result to standard normal RVs 
    MonteCarlo::ConvertLHStoStdNorm<size_t,Float> ( LHSResult );


    size_t maxP = 8;
    auto indices = BasisFunctions::MultiIndex ( dim, maxP );

    MatrixXF Correl ( 2, 2 );

    Correl ( 0, 0 ) = 1.0;
    Correl ( 1, 1 ) = 1.0;


    auto StdNormRVs = MonteCarlo::CombineRVs<Float,Complex> ( 

            Correl, 
            LHSResult, 
            dim 

    );



    // Create a vector of ICDF functions for inverse transform 
    std::vector< std::function < Float(Float) > > ICDFs;

    boost::math::lognormal dist1 ( Float(24.11948805), Float(0.099751) );

    ICDFs.emplace_back ( 

        [dist1](const auto m) { 
            return quantile ( dist1, m );
        }

    );

    boost::math::lognormal dist2 ( Float(7.822797571), Float(0.0499687) );

    ICDFs.emplace_back ( 

        [dist2](const auto m) { 
            return quantile ( dist2, m );
        }

    );

    auto RVs = MonteCarlo::GenerateRVs<size_t,Float> ( 

        LHSResult, 
        Correl, 
        ICDFs, 
        dim 

    );


    /**
      * Generate surrogate model
      * 1. Evaluate basis functions using generated random variables 
      * 2. Evaluate FRF model using generated random variables 
      * 3. Perform least square fitting 
      */

    // Evaluate basis functions 
    auto Basis = 
        BasisFunctions::HermitePolynomials<size_t, Float, Complex> (
            indices, StdNormRVs, dim 
        );


    // Evaluate FRF model and assemble into vector 
    auto model = []( const auto w, const auto& pars ) {
        return MonteCarlo::EvaluateFRF<Float, Complex> ( w, pars );
    };

    Float omega = 30;

    auto bStdVec = MonteCarlo::EvaluateModel <

        size_t,Float,Complex, decltype(model)

    > (

        omega, 
        RVs, 
        model, 
        dim 

    );

    auto nRow = bStdVec.size();
    auto nCol = indices.size() / dim;

    Eigen::Map<MatrixXC> A ( Basis.data(), nCol, nRow );

    Eigen::Map<VectorXC> b ( bStdVec.data(), bStdVec.size() );

    VectorXC x = ( A * A.transpose() ).ldlt().solve( A * b);

    std::cout << x << std::endl;


    return 0;
} // main 

