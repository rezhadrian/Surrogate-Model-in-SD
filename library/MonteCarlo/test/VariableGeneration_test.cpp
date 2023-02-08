/**
  * @file VariableGeneration_test.cpp
  *
  * @brief 
  * Tests of functions to generate correlated random variables 
  *
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "MonteCarlo.hpp" 
#include <gtest/gtest.h> 

TEST ( ConvertLHStoStdNorm, MockLHSResult ) {

    typedef double Float;

    std::vector<Float> LHSResult {

        0.5000, 0.5250, 0.5500, 0.5750, 0.6000, 0.6250, 0.6500, 0.6750, 
        0.7000, 0.7250, 0.7500, 0.7750, 0.8000, 0.8250, 0.8500, 0.8750, 
        0.5000, 0.4750, 0.4500, 0.4250, 0.4000, 0.3750, 0.3500, 0.3250,
        0.3000, 0.2750, 0.2500, 0.2250, 0.2000, 0.1750, 0.1500, 0.1250,
        0.9210, 0.9215, 0.9250, 0.9300, 0.9350, 0.9400, 0.9450, 0.9500,
        0.9550, 0.9600, 0.9650, 0.9700, 0.9750, 0.9800, 0.9850, 0.9900,
        0.9950, 0.9960, 0.9970, 0.9980, 0.9990, 0.9995, 0.9996, 0.9997,
        0.9998, 0.9999, 0.0790, 0.0785, 0.0750, 0.0700, 0.0650, 0.0600,
        0.0550, 0.0500, 0.0450, 0.0400, 0.0350, 0.0300, 0.0250, 0.0200,
        0.0150, 0.0100, 0.0050, 0.0040, 0.0030, 0.0020, 0.0010, 0.0005,
        0.0004, 0.0003, 0.0002, 0.0001

    };

    std::vector<Float> StdNorm {

         0.000000,  0.062707,  0.125661,  0.189118,  0.253347,  0.318639,
         0.385320,  0.453762,  0.524401,  0.597760,  0.674490,  0.755415,
         0.841621,  0.934589,  1.036433,  1.150349,  0.000000, -0.062707,
        -0.125661, -0.189118, -0.253347, -0.318639, -0.385320, -0.453762,
        -0.524401, -0.597760, -0.674490, -0.755415, -0.841621, -0.934589,
        -1.036433, -1.150349,  1.411830,  1.415234,  1.439531,  1.475791,
         1.514102,  1.554774,  1.598193,  1.644854,  1.695398,  1.750686,
         1.811911,  1.880794,  1.959964,  2.053749,  2.170090,  2.326348,
         2.575829,  2.652070,  2.747781,  2.878162,  3.090232,  3.290527,
         3.352795,  3.431614,  3.540084,  3.719016, -1.411830, -1.415234,
        -1.439531, -1.475791, -1.514102, -1.554774, -1.598193, -1.644854,
        -1.695398, -1.750686, -1.811911, -1.880794, -1.959964, -2.053749,
        -2.170090, -2.326348, -2.575829, -2.652070, -2.747781, -2.878162,
        -3.090232, -3.290527, -3.352795, -3.431614, -3.540084, -3.719016

    };

    MonteCarlo::ConvertLHStoStdNorm<size_t,Float> ( LHSResult );

    Float tol = 0.000001;

    for ( auto i = 0; i < LHSResult.size(); i++ ) {

        EXPECT_NEAR ( LHSResult[i], StdNorm[i], tol);
        
    }

}

TEST ( CombineRVs, NoCorrelation ) {

    typedef double Float;
    typedef std::complex<Float> Complex;
    typedef std::vector<Float> Vector; 
    typedef Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> MatrixXF;
    typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> MatrixXC;

    size_t dim = 3;

    Vector RVs {

        1.0, 2.0, 3.0, 
        4.0, 5.0, 6.0, 
        7.0, 8.0, 9.0 

    };

    MatrixXF Correl ( dim, dim );

    Correl ( 0, 0 ) = 1.0;
    Correl ( 1, 1 ) = 1.0;
    Correl ( 2, 2 ) = 1.0;

    std::vector<Complex> result = MonteCarlo::CombineRVs<Float,Complex> (

        Correl, 
        RVs, 
        dim

    );


    ASSERT_EQ ( result.size(), RVs.size() );

    for ( auto i = 0; i < result.size(); i++ ) {

        EXPECT_EQ ( result[i].real(), RVs[i] );

    }

}

TEST ( CombineRVs, SomeCorrelation ) {

    typedef double Float;
    typedef std::complex<Float> Complex;
    typedef std::vector<Float> Vector; 
    typedef Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> MatrixXF;
    typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> MatrixXC;

    size_t dim = 3;

    Vector RVs {

        1.0, 2.0, 3.0, 
        4.0, 5.0, 6.0, 
        7.0, 8.0, 9.0 

    };

    MatrixXF Correl ( dim, dim );

    Correl ( 0, 0 ) = 1.0;
    Correl ( 1, 0 ) = 0.3;
    Correl ( 1, 1 ) = 1.0;
    Correl ( 2, 2 ) = 1.0;

    std::vector<Complex> result = MonteCarlo::CombineRVs<Float,Complex> (

        Correl, 
        RVs, 
        dim

    );

    Vector expected {

        1.000000000000000, 2.207878402833891, 3.000000000000000, 
        4.000000000000000, 5.969696007084728, 6.000000000000000, 
        7.000000000000000, 9.731513611335565, 9.000000000000000 

    };


    ASSERT_EQ ( result.size(), expected.size() );

    Float tol = 1e-12;
    for ( auto i = 0; i < result.size(); i++ ) {

        EXPECT_NEAR ( result[i].real(), expected[i], tol );

    }

}

TEST ( GenerateRVs, NoCorrelationStdNorm ) {

    typedef double Float;
    typedef std::complex<Float> Complex;
    typedef std::vector<Float> Vector; 
    typedef Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> MatrixXF;
    typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> MatrixXC;

    size_t dim = 2;

    Vector RVs {

        1.0, 2.0, 
        1.0, 2.0, 
        1.0, 2.0

    };

    MatrixXF Correl ( dim, dim );

    Correl ( 0, 0 ) = 1.0;
    Correl ( 1, 1 ) = 1.0;

    boost::math::normal dist ( Float(0.0), Float(1.0) );

    std::vector< std::function < Float(Float) > > ICDFs;

    ICDFs.emplace_back ( 

        [dist](const auto m) { 
            return quantile ( dist, m );
        }

    );

    ICDFs.emplace_back ( 

        [dist](const auto m) { 
            return quantile ( dist, m );
        }

    );

    auto result = MonteCarlo::GenerateRVs<size_t,Float> (

        RVs, 
        Correl, 
        ICDFs, 
        dim

    );

    ASSERT_EQ ( RVs.size(), result.size() );

    Float tol = 1e-12;

    for ( auto i = 0; i < result.size(); i++ ) {

        EXPECT_NEAR ( RVs[i], result[i], tol );

    }

}

TEST ( GenerateRVs, SomeCorrelationStdNorm ) {

    typedef double Float;
    typedef std::complex<Float> Complex;
    typedef std::vector<Float> Vector; 
    typedef Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> MatrixXF;
    typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> MatrixXC;

    size_t dim = 2;

    Vector RVs {

        1.0, 2.0, 
        1.0, 2.0, 
        1.0, 2.0

    };

    MatrixXF Correl ( dim, dim );

    Correl ( 0, 0 ) = 1.0;
    Correl ( 1, 0 ) = 0.3;
    Correl ( 1, 1 ) = 1.0;

    boost::math::normal dist ( Float(0.0), Float(1.0) );

    std::vector< std::function < Float(Float) > > ICDFs;

    ICDFs.emplace_back ( 

        [dist](const auto m) { 
            return quantile ( dist, m );
        }

    );

    ICDFs.emplace_back ( 

        [dist](const auto m) { 
            return quantile ( dist, m );
        }

    );

    auto result = MonteCarlo::GenerateRVs<size_t,Float> (

        RVs, 
        Correl, 
        ICDFs, 
        dim

    );

    Vector expected {

        1.000000000000000, 2.207878402833891, 
        1.000000000000000, 2.207878402833891, 
        1.000000000000000, 2.207878402833891

    };

    ASSERT_EQ ( expected.size(), result.size() );

    Float tol = 1e-12;

    for ( auto i = 0; i < result.size(); i++ ) {

        EXPECT_NEAR ( expected[i], result[i], tol );

    }

}

