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

TEST ( InvLogNormCDF_test, NormalInput ) {

    double quantile = 0.3;
    double mu = 0.4;
    double sig = 0.4;

    double result = MonteCarlo::InvLogNormCDF ( quantile, mu, sig );

    double expected = 1.209539604333566;
    double tol      = 0.000000000000001;

    EXPECT_NEAR ( result, expected, tol );


}

