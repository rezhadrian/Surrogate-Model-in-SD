/**
  * @file MultiIndex_test.cpp
  *
  * @brief test output of MultiIndex function 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "BasisFunctions.hpp" 
#include <gtest/gtest.h> 

TEST ( MultiIndex, Constant1D ) {

    std::vector<size_t> expected { 0 };

    auto result = BasisFunctions::MultiIndex ( 1, 0 );

    ASSERT_EQ ( result.size(), expected.size() );
    for ( auto i = 0; i < expected.size(); i++ ) {
        EXPECT_EQ ( result[i],  expected[i] );
    }

}

TEST ( MultiIndex, Constant2D ) {

    std::vector<size_t> expected { 0, 0 };

    auto result = BasisFunctions::MultiIndex ( 2, 0 );

    ASSERT_EQ ( result.size(), expected.size() );
    for ( auto i = 0; i < expected.size(); i++ ) {
        EXPECT_EQ ( result[i],  expected[i] );
    }

}

TEST ( MultiIndex, Constant3D ) {

    std::vector<size_t> expected { 0, 0, 0 };

    auto result = BasisFunctions::MultiIndex ( 3, 0 );

    ASSERT_EQ ( result.size(), expected.size() );
    for ( auto i = 0; i < expected.size(); i++ ) {
        EXPECT_EQ ( result[i],  expected[i] );
    }

}

TEST ( MultiIndex, Linear1D ) {

    std::vector<size_t> expected { 0, 1 };

    auto result = BasisFunctions::MultiIndex ( 1, 1 );

    ASSERT_EQ ( result.size(), expected.size() );
    for ( auto i = 0; i < expected.size(); i++ ) {
        EXPECT_EQ ( result[i],  expected[i] );
    }

}

TEST ( MultiIndex, Linear2D ) {

    std::vector<size_t> expected { 0, 0, 1, 0, 0, 1 };

    auto result = BasisFunctions::MultiIndex ( 2, 1 );

    ASSERT_EQ ( result.size(), expected.size() );
    for ( auto i = 0; i < expected.size(); i++ ) {
        EXPECT_EQ ( result[i],  expected[i] );
    }

}

TEST ( MultiIndex, Linear3D ) {

    std::vector<size_t> expected {
        0, 0, 0, 
        1, 0, 0,
        0, 1, 0,
        0, 0, 1 
    };

    auto result = BasisFunctions::MultiIndex ( 3, 1 );

    ASSERT_EQ ( result.size(), expected.size() );
    for ( auto i = 0; i < expected.size(); i++ ) {
        EXPECT_EQ ( result[i],  expected[i] );
    }

}

TEST ( MultiIndex, Quadratic1D ) {

    std::vector<size_t> expected { 0, 1, 2 };

    auto result = BasisFunctions::MultiIndex ( 1, 2 );

    ASSERT_EQ ( result.size(), expected.size() );
    for ( auto i = 0; i < expected.size(); i++ ) {
        EXPECT_EQ ( result[i],  expected[i] );
    }

}

TEST ( MultiIndex, Quadratic2D ) {

    std::vector<size_t> expected {
        0, 0,
        1, 0,
        0, 1,
        2, 0,
        1, 1,
        0, 2
    };

    auto result = BasisFunctions::MultiIndex ( 2, 2 );

    ASSERT_EQ ( result.size(), expected.size() );
    for ( auto i = 0; i < expected.size(); i++ ) {
        EXPECT_EQ ( result[i],  expected[i] );
    }

}

TEST ( MultiIndex, Quadratic3D ) {

    std::vector<size_t> expected {
        0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,
        2, 0, 0,
        1, 1, 0,
        1, 0, 1,
        0, 2, 0,
        0, 1, 1,
        0, 0, 2
    };

    auto result = BasisFunctions::MultiIndex ( 3, 2 );

    ASSERT_EQ ( result.size(), expected.size() );
    for ( auto i = 0; i < expected.size(); i++ ) {
        EXPECT_EQ ( result[i],  expected[i] );
    }

}

TEST ( MultiIndex, NegativePower ) {

    EXPECT_THROW ({

        try {

            auto result = BasisFunctions::MultiIndex ( 3, -1 );

        } catch ( const std::exception& e ) {

            EXPECT_STREQ (
                "MultiIndex: max order cannot be negative",
                e.what()
            );

            throw;

        }


    }, std::runtime_error );

}

TEST ( MultiIndex, ZeroDimension ) {

    EXPECT_THROW ({

        try {

            auto result = BasisFunctions::MultiIndex ( 0, 5 );

        } catch ( const std::exception& e ) {

            EXPECT_STREQ (
                "MultiIndex: dimension must be positive integers",
                e.what()
            );

            throw;

        }


    }, std::runtime_error );

}

TEST ( MultiIndex, NegativeDimension ) {

    EXPECT_THROW ({

        try {

            auto result = BasisFunctions::MultiIndex ( -2, 5 );

        } catch ( const std::exception& e ) {

            EXPECT_STREQ (
                "MultiIndex: dimension must be positive integers",
                e.what()
            );

            throw;

        }


    }, std::runtime_error );

}

