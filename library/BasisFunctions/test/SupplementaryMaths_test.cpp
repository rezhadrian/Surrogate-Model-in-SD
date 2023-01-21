/**
  * This file is part of a Surrogate Model library for structural dynamics. 
  *
  * The library is free: you can distribute it and/or modify it.
  * The library is distributed in the hope that it can be useful and helpful,
  * particularly for students and self - learning individuals.
  * 
  * The library comes WITHOUT ANY WARRANTY: without even the implied warranty 
  * of MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE.
  *
  * The library is licensed under the GNU General Public License v3.0
  *
  *
  * The library is based on works by:
  * - Felix Schneider
  * - Iason Papaioannou
  * 
  * Chair of Structural Mechanics 
  * Engineering Risk Analysis Group 
  *
  * Technische Universitaet Muenchen
  *
  * www.cee.ed.tum.de/bm
  * www.cee.ed.tum.de/era 
  */ 

/**
  * @file SupplementaryMaths_test.cpp
  *
  * @brief test supplementary maths functions implementations 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "SupplementaryMaths.hpp"  
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

