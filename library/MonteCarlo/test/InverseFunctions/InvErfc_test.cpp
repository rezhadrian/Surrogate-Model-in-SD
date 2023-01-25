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
        0.6, 0.7, 0.8, 0.9

    };

    std::vector <double> xL( p.size() );
    std::vector <double> xR ( p.size() );

    for ( auto i = 0; i < p.size(); i++ ) {

        xL[i] = MonteCarlo::InvErfc (       p[i] );
        xR[i] = MonteCarlo::InvErfc ( 2.0 - p[i] );

    }

    std::vector <double> expected {

        1.163087153676674, 0.906193802436823, 0.732869077959217, 
        0.595116081449995, 0.476936276204470, 0.370807158593558, 
        0.272462714726754, 0.179143454621292, 0.088855990494258, 
         
    };


    double tol = 0.000000000000001;

    for ( auto i = 0; i < p.size(); i++ ) {

        EXPECT_NEAR ( xL[i],  expected[i], tol);
        EXPECT_NEAR ( xR[i], -expected[i], tol);

    }

}

