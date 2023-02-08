/**
  * @file TripleHermite_test.cpp
  *
  * @brief test function 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "BasisFunctions.hpp" 
#include <gtest/gtest.h> 

TEST ( ETripleHermite, TypicalInput ) {

    typedef double Float;

    std::vector<size_t> input1 {

        0, 0, 0,
        0, 0, 0,
        0, 0, 0

    };

    std::vector<size_t> input2 {

        0, 0, 0,
        1, 1, 1, 
        2, 2, 2 

    };

    std::vector<size_t> input3 {

        0, 1, 2, 
        0, 1, 2, 
        0, 1, 2

    };

    std::vector<Float> expected {

        1.0, 0.0, 0.0, 
        0.0, 1.0, 0.0, 
        0.0, 0.0, 1.0 

    };

    double tol = 1e-12;

    for ( auto i = 0; i < expected.size(); i++ ) {

        auto result = BasisFunctions::ETripleHermite<size_t,Float> (
            input1[i], input2[i], input3[i]
        );

        EXPECT_NEAR ( result, expected[i], tol );

    }

}

TEST ( EMultiTripleHermite, TypicalInput ) {

    typedef double Float; 

    size_t dim = 3; 
    size_t k   = 0; 

    std::vector<size_t> indices {

        0, 1, 2, 
        5, 4, 3 

    };

    auto result = BasisFunctions::EMultiTripleHermite<size_t,Float> (

            indices, dim, k 

    );

    std::vector<Float> expected {

        1.0, 0.0, 
        0.0, 1.0 

    };

    ASSERT_EQ ( result.size(), expected.size() );

    Float tol = 1e-12;

    for ( auto i = 0; i < result.size(); i++ ) {

        EXPECT_NEAR ( result[i], expected[i], tol );

    }

}

