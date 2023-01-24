/**
  * @file InvErfc.cpp
  *
  * @brief test the function inverse for complement of error function 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "MonteCarlo.hpp" 
#include <gtest/gtest.h>

TEST ( InvErfc, ZeroToTwo ) {

    std::vector <double> p {

        0.1, 0.2, 0.3, 0.4, 0.5,
        0.6, 0.7, 0.8, 0.9, 1.0,
        1.1, 1.2, 1.3, 1.4, 1.5,
        1.6, 1.7, 1.8, 1.9

    };

    std::vector <double> x ( p.size() );

    for ( auto i = 0; i < p.size(); i++ ) {

        x[i] = MonteCarlo::InvErfc ( p[i] );

    }

    std::vector <double> expected {

         1.16308715,  0.90619380,  0.73286908,  0.59511608,
         0.47693628,  0.37080716,  0.27246272,  0.17914346,
         0.08885599,  0.00000000, -0.08885599, -0.17914346,
        -0.27246272, -0.37080716, -0.47693628, -0.59511608,
        -0.73286908, -0.90619380, -1.16308715

    };


    double tol = 0.00000001;

    for ( auto i = 0; i < p.size(); i++ ) {

        EXPECT_NEAR ( x[i], expected[i], tol);

    }

}

