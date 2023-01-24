/**
  * @file SupplementaryMaths_test.cpp
  *
  * @brief test supplementary maths functions implementations 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "BasisFunctions.hpp" 
#include <gtest/gtest.h> 

TEST ( BinomialCoefficient, TypicalInput ) {

    EXPECT_EQ ( BasisFunctions::Binomial<size_t> ( 4, 0 ), 1 );
    EXPECT_EQ ( BasisFunctions::Binomial<size_t> ( 4, 1 ), 4 );
    EXPECT_EQ ( BasisFunctions::Binomial<size_t> ( 4, 2 ), 6 );
    EXPECT_EQ ( BasisFunctions::Binomial<size_t> ( 4, 3 ), 4 );
    EXPECT_EQ ( BasisFunctions::Binomial<size_t> ( 4, 4 ), 1 );

    EXPECT_EQ ( BasisFunctions::Binomial<size_t> ( 3, 0 ), 1 );
    EXPECT_EQ ( BasisFunctions::Binomial<size_t> ( 3, 1 ), 3 );
    EXPECT_EQ ( BasisFunctions::Binomial<size_t> ( 3, 2 ), 3 );
    EXPECT_EQ ( BasisFunctions::Binomial<size_t> ( 3, 3 ), 1 );

}

TEST ( BinomialCoefficient, BadInput ) {

    EXPECT_EQ ( BasisFunctions::Binomial<size_t> (  7,  9 ), 0 );
    EXPECT_EQ ( BasisFunctions::Binomial<size_t> (  0,  0 ), 1 );
    EXPECT_EQ ( BasisFunctions::Binomial<size_t> (  0,  2 ), 0 );
    EXPECT_EQ ( BasisFunctions::Binomial<size_t> (  9, 12 ), 0 );
    EXPECT_EQ ( BasisFunctions::Binomial<size_t> (  4,  5 ), 0 );

}

TEST ( BinomialCoefficient, UnacceptableInput ) {

    EXPECT_THROW ({

        try {
            
            auto output = BasisFunctions::Binomial<int> (-1,0);

        } catch ( const std::exception& e ) {

            EXPECT_STREQ (
                "Binomial: currently doesn't support negative values",
                e.what()
            );

            throw;

        }
            
    }, std::runtime_error );

    EXPECT_THROW ({

        try {
            
            auto output = BasisFunctions::Binomial<int> (0,-1);

        } catch ( const std::exception& e ) {

            EXPECT_STREQ (
                "Binomial: currently doesn't support negative values",
                e.what()
            );

            throw;

        }
            
    }, std::runtime_error );

}

