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
TEST ( EvaluateModel, StandardEandMu ) {

    typedef std::complex<long double> Complex;
    typedef std::vector<long double> Vector;

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
         decltype(model)> ( omega, params, model, 2 );

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

