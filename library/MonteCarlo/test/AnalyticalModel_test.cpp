/**
  * @file AnalyticalModel_test.cpp
  *
  * @brief test analytical Frequency Response Function (FRF) models 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "MonteCarlo.hpp" 
#include <gtest/gtest.h> 

TEST ( EvaluateFRF, StandardEandMu ) {

    typedef std::complex<long double> Complex;
    typedef std::vector<long double> Vector;

    Vector params {

        3e10, 2500 

    };

    long double omega = 10;

    Complex FRF;

    FRF = MonteCarlo::EvaluateFRF<long double,Complex> ( omega, params );


    Complex expected ( 0.091491467231740, -0.001156587941292 );

    long double tol =  0.000000000000001;

    EXPECT_NEAR ( FRF.real(), expected.real(), tol );
    EXPECT_NEAR ( FRF.imag(), expected.imag(), tol );

}

TEST ( EvaluateFRF, NegativeE ) {

    typedef std::complex<long double> Complex;
    typedef std::vector<long double> Vector;

    Vector params {

        -3e10, 2500 

    };

    long double omega = 10;

    Complex FRF;

    try {

        FRF = MonteCarlo::EvaluateFRF<long double,Complex> ( omega, params );

    } catch ( const std::exception& e ) {

        EXPECT_STREQ (
            "EvaluateFRF: Elasticity modulus must be positive", 
            e.what()
        );

    }

}

TEST ( EvaluateFRF, NegativeMu ) {

    typedef std::complex<long double> Complex;
    typedef std::vector<long double> Vector;

    Vector params {

        3e10, -2500 

    };

    long double omega = 10;

    Complex FRF;

    try {

        FRF = MonteCarlo::EvaluateFRF<long double,Complex> ( omega, params );

    } catch ( const std::exception& e ) {

        EXPECT_STREQ (
            "EvaluateFRF: Mu must be positive", 
            e.what()
        );

    }

}

TEST ( EvaluateModel, StandardEandMu ) {

    typedef std::complex<long double> Complex;
    typedef std::vector<long double> Vector;

    size_t dim = 2;

    Vector params {

        3.0e10, 2500, 
        2.0e10, 3000, 
        2.5e10, 2800

    };

    long double omega = 10;


    auto model = []( const auto w, const auto& pars ) {
        return MonteCarlo::EvaluateFRF<long double, Complex> ( w, pars );
    };


    auto result = MonteCarlo::EvaluateModel<size_t, long double, Complex, 
         decltype(model)> ( omega, params, model, dim );

    std::vector<Complex> expected {

        Complex ( 0.091491467231740, -0.001156587941292 ),
        Complex ( 0.177661504804962, -0.003251201702618 ), 
        Complex ( 0.126952244180099, -0.001921005625179 )

    };

    ASSERT_EQ ( result.size(), expected.size() );

    long double tol =  0.000000000000001;

    for ( auto i = 0; i < result.size(); i++ ) {

        EXPECT_NEAR ( result[i].real(), expected[i].real(), tol );
        EXPECT_NEAR ( result[i].imag(), expected[i].imag(), tol );

    }

}

TEST ( EvaluateModel, WrongDim ) {

    typedef std::complex<long double> Complex;
    typedef std::vector<long double> Vector;

    size_t dim = 0;

    Vector params {

        3.0e10, 2500, 
        2.0e10, 3000, 
        2.5e10, 2800

    };

    long double omega = 10;


    auto model = []( const auto w, const auto& pars ) {
        return MonteCarlo::EvaluateFRF<long double, Complex> ( w, pars );
    };


    try {

        auto result = MonteCarlo::EvaluateModel<size_t, long double, 
             Complex, decltype(model)> ( omega, params, model, dim );

    } catch ( const std::exception& e ) {

        EXPECT_STREQ (
            "EvaluateModel: dimension must be positive", 
            e.what()
        );

    }

}

