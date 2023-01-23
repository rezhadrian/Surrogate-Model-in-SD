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
  * @file TotalTruncations_test.cpp
  *
  * @brief test functions TotalTruncation 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "Truncations.hpp" 
#include <gtest/gtest.h> 

TEST ( TotalTruncations, NoTruncations ) {

    std::vector <size_t> indices {

        0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        0, 0, 1

    };

    std::vector <size_t> expected {

        0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        0, 0, 1

    };

    size_t  dim = 3;
    size_t pMax = 1;

    BasisFunctions::TotalTruncation ( indices, dim, pMax );

    ASSERT_EQ ( indices.size(), expected.size() );

    for ( auto i = 0; i < expected.size(); i++ ) {

        EXPECT_EQ ( indices[i], expected[i] );

    }

}

TEST ( TotalTruncations, ZeroPMax ) {

    std::vector <size_t> indices {

        0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        0, 0, 1

    };

    std::vector <size_t> expected {

        0, 0, 0

    };

    size_t  dim = 3;
    size_t pMax = 0;

    BasisFunctions::TotalTruncation ( indices, dim, pMax );

    ASSERT_EQ ( indices.size(), expected.size() );

    for ( auto i = 0; i < expected.size(); i++ ) {

        EXPECT_EQ ( indices[i], expected[i] );

    }

}

TEST ( TotalTruncations, FinitePMax ) {

    std::vector <size_t> indices {

        0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,
        2, 0, 0,
        1, 1, 0,
        1, 0, 1,
        0, 2, 0

    };

    std::vector <size_t> expected {

        0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,

    };

    size_t  dim = 3;
    size_t pMax = 1;

    BasisFunctions::TotalTruncation ( indices, dim, pMax );

    ASSERT_EQ ( indices.size(), expected.size() );

    for ( auto i = 0; i < expected.size(); i++ ) {

        EXPECT_EQ ( indices[i], expected[i] );

    }

}

TEST ( TotalTruncations, WrongDimension ) {

    std::vector <size_t> indices {

        0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,
        2, 0, 0,
        1, 1, 0,
        1, 0, 1,
        0, 2, 0

    };

    size_t  dim = 7;
    size_t pMax = 1;

    EXPECT_THROW ({

        try {
            
            BasisFunctions::TotalTruncation ( indices, dim, pMax );

        } catch ( const std::exception& e ) {

            EXPECT_STREQ (
                "Total truncation: Indices size not a multiple of dimension",
                e.what()
            );

            throw;

        }
            
    }, std::runtime_error );

}

