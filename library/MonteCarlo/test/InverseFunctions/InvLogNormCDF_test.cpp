/**
  * @file InvLogNormCDF_test.cpp 
  *
  * @brief test inverse function of log normal CDF 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "MonteCarlo.hpp" 
#include <gtest/gtest.h> 

TEST ( InvLogNormCDF_test, NormalInput1 ) {

    double mu  = 0.0;
    double sig = 1.0;

    std::vector<double> quantiles {

        0.1, 0.2, 0.3, 0.4, 0.5, 
        0.6, 0.7, 0.8, 0.9

    };

    std::vector<double> results ( quantiles.size(), 0.0 );

    for ( auto i = 0; i < results.size(); i++ ) {

        results[i] = MonteCarlo::InvLogNormCDF ( quantiles[i], mu, sig );

    }

    std::vector<double> expected {

        0.277606241852010, 0.431011186881839, 0.591910100609554, 
        0.776198414156351, 1.000000000000000, 1.288330382750007, 
        1.689445743484004, 2.320125394504318, 3.602224479279158 

    };

    double tol = 0.00000000000005;

    for ( auto i = 0; i < results.size(); i++ ) {

        EXPECT_NEAR ( results[i], expected[i], tol );

    }

}

TEST ( InvLogNormCDF_test, NormalInput2 ) {

    double mu  = 1.0;
    double sig = 0.1;

    std::vector<double> quantiles {

        0.1, 0.2, 0.3, 0.4, 0.5, 
        0.6, 0.7, 0.8, 0.9

    };

    std::vector<double> results ( quantiles.size(), 0.0 );

    for ( auto i = 0; i < results.size(); i++ ) {

        results[i] = MonteCarlo::InvLogNormCDF ( quantiles[i], mu, sig );

    }

    std::vector<double> expected {

        2.391318394729171, 2.498868118230018, 2.579408086382776, 
        2.650279986463947, 2.718281828459045, 2.788028486299391, 
        2.864632447242061, 2.956961211768318, 3.089950763234730 

    };

    double tol = 0.00000000000005;

    for ( auto i = 0; i < results.size(); i++ ) {

        EXPECT_NEAR ( results[i], expected[i], tol );

    }

}

