/**
  * @file MassSpringDamper_test.cpp 
  *
  * @brief 
  * Tests of the Mass-Spring-Damper analytical SD model 
  * 
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "AnalyticalModel.hpp" 
#include <gtest/gtest.h> 

TEST ( MassSpringDamper, NormalModel ) {

    typedef double Float;
    typedef std::vector<Float> Vector;
    typedef std::complex<Float> Complex; 

    typedef std::vector<Complex> VectorC;

    typedef Analytical::MassSpringDamper<size_t,Float,Complex> Model;
    // typedef Analytical::Model<size_t,Float,Complex> AbsModel;

    Vector Masses  { 2.0, 1.0 };
    Vector Dampers { 0.0, 0.0 };
    Vector Springs { 1.0, 1.0 };

    Vector AdditionalSprings { 1.0, -2.0 };

    Float omega = 0.0;

    Model SDModel ( Masses, Dampers, Springs );


    Vector MM {

        2.0, 0.0,
        0.0, 1.0 

    };

    Vector CM {

        0.0, 0.0, 
        0.0, 0.0 

    };

    Vector KM {

         2.0, -1.0, 
        -1.0,  1.0

    };

    Vector AKM {

        -1.0,  2.0, 
         2.0, -2.0 

    };

    VectorC Force { 1.0, 1.0 };

    auto start = AdditionalSprings.begin();
    auto end   = AdditionalSprings.end();
    auto result = SDModel.ComputeResponse ( 
        Force, omega, start, end );

    VectorC expected { 1.0, 0.0 };

    // auto DS = SDModel.DynamicStiffness ( omega, AdditionalSprings );


    Float tol = 1e-12;

    ASSERT_EQ ( result.size(), expected.size() );
    for ( auto i = 0; i < expected.size(); i++ ) {

        EXPECT_NEAR ( result[i].real(), expected[i].real(), tol );
        EXPECT_NEAR ( result[i].imag(), expected[i].imag(), tol );

    }

}

