/**
  * @file LatinHypercubeSampling_test.cpp
  *
  * @brief test functions to produce latin hypercube sample
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "MonteCarlo.hpp" 
#include <gtest/gtest.h> 

TEST ( LatinHypercubeSampling, SingleIntervalSingleSample ) {

    size_t nInterval = 1;
    size_t nSample   = 1;

    auto result = 
        MonteCarlo::LHS <

            size_t, double,
            std::random_device,
            std::default_random_engine,
            std::uniform_real_distribution<double>,
            std::mt19937

        > ( nInterval, nSample );

    ASSERT_EQ ( result.size(), nInterval * nSample );

    for ( auto i = 0; i < result.size(); i++ ) {

        EXPECT_TRUE ( result[i] >= 0.0 && result[i] <= 1.0 );

    }

}

TEST ( LatinHypercubeSampling, MultipleIntervalsSingleSample ) {

    size_t nInterval = 40;
    size_t nSample   = 1;

    auto result = 
        MonteCarlo::LHS <

            size_t, double,
            std::random_device,
            std::default_random_engine,
            std::uniform_real_distribution<double>,
            std::mt19937

        > ( nInterval, nSample );

    double range = 1.0 / nInterval;

    ASSERT_EQ ( result.size(), nInterval * nSample );

    for ( auto i = 0; i < nSample; i++ ) {

        std::vector <bool> indices ( nInterval, false );

        for ( auto j = 0; j < nInterval; j++ ) {
            
            for ( auto k = 0; k < nInterval; k++ ) {

                if (
                    result[i*nInterval + j] >=  k    * range &&
                    result[i*nInterval + j] <= (k+1) * range
                ) {
                    indices[k] = true;
                }

            }

        }

        EXPECT_TRUE ( 

            std::all_of ( 

                indices.begin(),
                indices.end(),

                []( auto m ) { return m; }
                        
            ) 

        );

    }

}

TEST ( LatinHypercubeSampling, SingleIntervalMultipleSamples ) {

    size_t nInterval = 1;
    size_t nSample   = 50;

    auto result = 
        MonteCarlo::LHS <

            size_t, double,
            std::random_device,
            std::default_random_engine,
            std::uniform_real_distribution<double>,
            std::mt19937

        > ( nInterval, nSample );

    ASSERT_EQ ( result.size(), nInterval * nSample );

    for ( auto i = 0; i < result.size(); i++ ) {

        EXPECT_TRUE ( result[i] >= 0.0 && result[i] <= 1.0 );

    }

}

TEST ( LatinHypercubeSampling, MultipleIntervalsMultipleSamples ) {

    size_t nInterval = 40;
    size_t nSample   = 50;

    auto result = 
        MonteCarlo::LHS <

            size_t, double,
            std::random_device,
            std::default_random_engine,
            std::uniform_real_distribution<double>,
            std::mt19937

        > ( nInterval, nSample );

    double range = 1.0 / nInterval;

    ASSERT_EQ ( result.size(), nInterval * nSample );

    for ( auto i = 0; i < nSample; i++ ) {

        std::vector <bool> indices ( nInterval, false );

        for ( auto j = 0; j < nInterval; j++ ) {
            
            for ( auto k = 0; k < nInterval; k++ ) {

                if (
                    result[i*nInterval + j] >=  k    * range &&
                    result[i*nInterval + j] <= (k+1) * range
                ) {
                    indices[k] = true;
                }

            }

        }

        EXPECT_TRUE ( 

            std::all_of ( 

                indices.begin(),
                indices.end(),

                []( auto m ) { return m; }
                        
            ) 

        );

    }

}

