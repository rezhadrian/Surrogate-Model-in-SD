/**
  * @file InverseCDFs_test.cpp
  *
  * @brief test functions to compute inverse CDF of RV vector 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "MonteCarlo.hpp" 
#include <gtest/gtest.h> 

TEST ( InverseCDFs, StdNormRightCentral ) {

    std::vector<double> quantiles {

        0.5000, 0.5250,
        0.5500, 0.5750,
        0.6000, 0.6250,
        0.6500, 0.6750,
        0.7000, 0.7250,
        0.7500, 0.7750,
        0.8000, 0.8250,
        0.8500, 0.8750,
        0.9000

    };

    std::vector<double> expected {

        0.000000, 0.062707,
        0.125661, 0.189118,
        0.253347, 0.318639,
        0.385320, 0.453762,
        0.524401, 0.597760,
        0.674490, 0.755415,
        0.841621, 0.934589,
        1.036433, 1.150349,
        1.281552

    };

    double tol = 0.000001;

    // choose inverse CDF : inverse standard normal CDF
    auto ICDF = [](const auto m){
        return MonteCarlo::InvStdNormCDF<size_t,double> ( m ) ;
    };

    // apply ICDF to all RVs in the vector 
    MonteCarlo::InverseCDFs<size_t,double,decltype(ICDF)>(quantiles,ICDF);

    for ( auto i = 0; i < quantiles.size(); i++ ) {

        EXPECT_NEAR ( quantiles[i], expected[i], tol);
        
    }

}

TEST ( InverseCDFs, StdNormLeftCentral ) {

    std::vector<double> quantiles {

        0.5000, 0.4750,
        0.4500, 0.4250,
        0.4000, 0.3750,
        0.3500, 0.3250,
        0.3000, 0.2750,
        0.2500, 0.2250,
        0.2000, 0.1750,
        0.1500, 0.1250,
        0.1000

    };

    std::vector<double> expected {

         0.000000, -0.062707,
        -0.125661, -0.189118,
        -0.253347, -0.318639,
        -0.385320, -0.453762,
        -0.524401, -0.597760,
        -0.674490, -0.755415,
        -0.841621, -0.934589,
        -1.036433, -1.150349,
        -1.281552

    };

    double tol = 0.000001;

    // choose inverse CDF : inverse standard normal CDF
    auto ICDF = [](const auto m){
        return MonteCarlo::InvStdNormCDF<size_t,double> ( m ) ;
    };

    // apply ICDF to all RVs in the vector 
    MonteCarlo::InverseCDFs<size_t,double,decltype(ICDF)>(quantiles,ICDF);

    for ( auto i = 0; i < quantiles.size(); i++ ) {

        EXPECT_NEAR ( quantiles[i], expected[i], tol);
        
    }

}

TEST ( InverseCDFs, StdNormRightTail ) {

    std::vector<double> quantiles {

        0.9210, 0.9215,
        0.9250, 0.9300,
        0.9350, 0.9400,
        0.9450, 0.9500,
        0.9550, 0.9600,
        0.9650, 0.9700,
        0.9750, 0.9800,
        0.9850, 0.9900,
        0.9950, 0.9960,
        0.9970, 0.9980,
        0.9990, 0.9995,
        0.9996, 0.9997,
        0.9998, 0.9999

    };

    std::vector<double> expected {

        1.411830, 1.415234,
        1.439531, 1.475791,
        1.514102, 1.554774,
        1.598193, 1.644854,
        1.695398, 1.750686,
        1.811911, 1.880794,
        1.959964, 2.053749,
        2.170090, 2.326348,
        2.575829, 2.652070,
        2.747781, 2.878162,
        3.090232, 3.290527,
        3.352795, 3.431614,
        3.540084, 3.719016

    };

    double tol = 0.000001;

    // choose inverse CDF : inverse standard normal CDF
    auto ICDF = [](const auto m){
        return MonteCarlo::InvStdNormCDF<size_t,double> ( m ) ;
    };

    // apply ICDF to all RVs in the vector 
    MonteCarlo::InverseCDFs<size_t,double,decltype(ICDF)>(quantiles,ICDF);

    for ( auto i = 0; i < quantiles.size(); i++ ) {

        EXPECT_NEAR ( quantiles[i], expected[i], tol);
        
    }

}

TEST ( InverseCDFs, StdNormLeftTail ) {

    std::vector<double> quantiles {

        0.0790, 0.0785,
        0.0750, 0.0700,
        0.0650, 0.0600,
        0.0550, 0.0500,
        0.0450, 0.0400,
        0.0350, 0.0300,
        0.0250, 0.0200,
        0.0150, 0.0100,
        0.0050, 0.0040,
        0.0030, 0.0020,
        0.0010, 0.0005,
        0.0004, 0.0003,
        0.0002, 0.0001

    };

    std::vector<double> expected {

        -1.411830, -1.415234,
        -1.439531, -1.475791,
        -1.514102, -1.554774,
        -1.598193, -1.644854,
        -1.695398, -1.750686,
        -1.811911, -1.880794,
        -1.959964, -2.053749,
        -2.170090, -2.326348,
        -2.575829, -2.652070,
        -2.747781, -2.878162,
        -3.090232, -3.290527,
        -3.352795, -3.431614,
        -3.540084, -3.719016

    };

    double tol = 0.000001;

    // choose inverse CDF : inverse standard normal CDF
    auto ICDF = [](const auto m){
        return MonteCarlo::InvStdNormCDF<size_t,double> ( m ) ;
    };

    // apply ICDF to all RVs in the vector 
    MonteCarlo::InverseCDFs<size_t,double,decltype(ICDF)>(quantiles,ICDF);

    for ( auto i = 0; i < quantiles.size(); i++ ) {

        EXPECT_NEAR ( quantiles[i], expected[i], tol);
        
    }

}

TEST ( ComputeCDFs, PositiveInput ) {

    std::vector<long double> X {

        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 
        0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 
        2.0, 2.1, 2.2, 2.3, 2.4, 2.5 

    };

    std::vector<long double> expected {

        0.500000000000000, 0.539827837277029, 0.579259709439103, 
        0.617911422188953, 0.655421741610324, 0.691462461274013, 
        0.758036347776927, 0.815939874653240, 0.864333939053617, 
        0.903199515414390, 0.933192798731142, 0.955434537241457, 
        0.977249868051821, 0.982135579437183, 0.986096552486501, 
        0.989275889978324, 0.991802464075404, 0.993790334674224
        
    };

    auto CDF = [](const auto m){
        return MonteCarlo::StdNormCDF<long double> ( m ) ;
    };

    // apply ICDF to all RVs in the vector 
    MonteCarlo::ComputeCDFs<size_t,long double,decltype(CDF)>(X,CDF);

    long double tol = 0.000000000000001; 

    for ( auto i = 0; i < X.size(); i++ ) {

        EXPECT_NEAR ( X[i], expected[i], tol);
        
    }

}

