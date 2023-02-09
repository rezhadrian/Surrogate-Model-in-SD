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

    typedef Analytical::MassSpringDamper<size_t,Float,Complex> Model;

    Vector Masses  { 2.0, 1.0 };
    Vector Dampers { 0.0, 0.0 };
    Vector Springs { 1.0, 1.0 };

    Model SDModel ( Masses, Dampers, Springs );

    auto MMatrix = SDModel.ComputeMassMatrix ();
    auto CMatrix = SDModel.ComputeDampingMatrix ();
    auto KMatrix = SDModel.ComputeStiffnessMatrix ();

    Vector MMExpected {

        2.0, 0.0,
        0.0, 1.0 

    };

    Vector CMExpected {

        0.0, 0.0, 
        0.0, 0.0 

    };

    Vector KMExpected {

         2.0, -1.0, 
        -1.0,  1.0

    };

    ASSERT_EQ ( MMatrix.size(), MMExpected.size() );
    ASSERT_EQ ( CMatrix.size(), CMExpected.size() );
    ASSERT_EQ ( KMatrix.size(), KMExpected.size() );

    Float tol = 1e-12;

    for ( auto i = 0; i < MMExpected.size(); i++ ) {

        EXPECT_NEAR ( MMatrix[i], MMExpected[i], tol );
        EXPECT_NEAR ( CMatrix[i], CMExpected[i], tol );
        EXPECT_NEAR ( KMatrix[i], KMExpected[i], tol );

    }

}

