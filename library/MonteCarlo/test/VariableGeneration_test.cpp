/**
  * @file VariableGeneration_test.cpp
  *
  * @brief test function to generate variable 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "MonteCarlo.hpp" 
#include <gtest/gtest.h> 

TEST ( CombineRVs, NoChange ) {

    typedef MonteCarlo::DenseLTriangularMatrix<size_t, double> Matrix;
    typedef std::vector<double> Vector;

    size_t dim = 6;

    Matrix A ( dim );

    for ( auto i = 0; i < dim; i++ ) {

        A(i,i) = 1.0;

    }

    Vector RVs {

        1.0, 7.0,
        2.0, 8.0, 
        3.0, 9.0, 
        4.0, 8.0, 
        5.0, 7.0, 
        6.0, 6.0 

    };

    Vector expected {

        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 
        7.0, 8.0, 9.0, 8.0, 7.0, 6.0 

    };

    Vector result = MonteCarlo::CombineRVs ( A, RVs, dim );

    ASSERT_EQ ( result.size(), expected.size() );

    for ( auto i = 0; i < result.size(); i++ ) {

        EXPECT_EQ ( result[i], expected[i] );

    }

}

TEST ( CombineRVs, SomeChanges ) {

    typedef MonteCarlo::DenseLTriangularMatrix<size_t, double> Matrix;
    typedef std::vector<double> Vector;

    size_t dim = 6;

    Matrix A ( dim );

    for ( auto i = 0; i < dim; i++ ) {

        A(i,i) = 1.0;

    }

    A(4,3) = 0.5;
    A(5,0) = 1.0;

    Vector RVs {

        1.0, 7.0,
        2.0, 8.0, 
        3.0, 9.0, 
        4.0, 8.0, 
        5.0, 7.0, 
        6.0, 6.0 

    };

    Vector expected {

        1.0, 2.0, 3.0, 4.0, 7.0, 7.0, 
        7.0, 8.0, 9.0, 8.0, 11.0, 13.0 

    };

    Vector result = MonteCarlo::CombineRVs ( A, RVs, dim );

    ASSERT_EQ ( result.size(), expected.size() );

    for ( auto i = 0; i < result.size(); i++ ) {

        EXPECT_EQ ( result[i], expected[i] );

    }

}

