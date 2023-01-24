/**
  * @file TotalTruncations_test.cpp
  *
  * @brief test functions TotalTruncation 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "BasisFunctions.hpp" 
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

