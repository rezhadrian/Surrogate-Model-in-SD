/**
  * @file LatinHypercubeSampling_test.cpp
  *
  * @brief 
  * Tests of functions to produce latin hypercube samples 
  *
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "MonteCarlo.hpp" 
#include <gtest/gtest.h> 

TEST ( UnsortedLHS, SingleIntervalSingleSample ) {

    size_t nPoints = 1;
    size_t dim     = 1;

    auto result = 
        MonteCarlo::UnsortedLHS < size_t, double > ( nPoints, dim );

    double range = 1.0 / nPoints;

    ASSERT_EQ ( result.size(), nPoints * dim );

    for ( auto j = 0; j < dim; j++ ) {

        std::vector <bool> indices ( nPoints, false );

        for ( auto i = 0; i < nPoints; i++ ) {
            
            for ( auto k = 0; k < nPoints; k++ ) {

                if (
                    result[i + j * nPoints] >=  k    * range &&
                    result[i + j * nPoints] <= (k+1) * range
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

TEST ( UnsortedLHS, MultipleIntervalsSingleSample ) {

    size_t nPoints = 40;
    size_t dim     = 1;

    auto result = 
        MonteCarlo::UnsortedLHS < size_t, double > ( nPoints, dim );

    double range = 1.0 / nPoints;

    ASSERT_EQ ( result.size(), nPoints * dim );

    for ( auto j = 0; j < dim; j++ ) {

        std::vector <bool> indices ( nPoints, false );

        for ( auto i = 0; i < nPoints; i++ ) {
            
            for ( auto k = 0; k < nPoints; k++ ) {

                if (
                    result[i + j * nPoints] >=  k    * range &&
                    result[i + j * nPoints] <= (k+1) * range
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

TEST ( UnsortedLHS, SingleIntervalsMultipleSample ) {

    size_t nPoints = 1;
    size_t dim     = 200;

    auto result = 
        MonteCarlo::UnsortedLHS < size_t, double > ( nPoints, dim );

    double range = 1.0 / nPoints;

    ASSERT_EQ ( result.size(), nPoints * dim );

    for ( auto j = 0; j < dim; j++ ) {

        std::vector <bool> indices ( nPoints, false );

        for ( auto i = 0; i < nPoints; i++ ) {
            
            for ( auto k = 0; k < nPoints; k++ ) {

                if (
                    result[i + j * nPoints] >=  k    * range &&
                    result[i + j * nPoints] <= (k+1) * range
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

TEST ( UnsortedLHS, MultipleIntervalsMultipleSamples ) {

    size_t nPoints = 40;
    size_t dim     = 50;

    auto result = 
        MonteCarlo::UnsortedLHS < size_t, double > ( nPoints, dim );

    double range = 1.0 / nPoints;

    ASSERT_EQ ( result.size(), nPoints * dim );

    for ( auto j = 0; j < dim; j++ ) {

        std::vector <bool> indices ( nPoints, false );

        for ( auto i = 0; i < nPoints; i++ ) {
            
            for ( auto k = 0; k < nPoints; k++ ) {

                if (
                    result[i + j * nPoints] >=  k    * range &&
                    result[i + j * nPoints] <= (k+1) * range
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

TEST ( UnsortedLHS, WrongNPoints ) {

    size_t nPoints = 0;
    size_t dim     = 2;

    try {

        auto result = MonteCarlo::UnsortedLHS<size_t,double> ( nPoints, dim );

    } catch ( const std::exception& e ) {

        EXPECT_STREQ (
            "LHS: Number of points must be positive integer",
            e.what()
        );

    }

}

TEST ( UnsortedLHS, WrongDim ) {


    size_t nPoints = 3;
    size_t dim     = 0;

    try {

        auto result = MonteCarlo::UnsortedLHS<size_t,double> ( nPoints, dim );

    } catch ( const std::exception& e ) {

        EXPECT_STREQ (
            "LHS: dimension must be positive integer",
            e.what()
        );

    }

}

TEST ( TransposerIndices, SingleIntervalSingleSample ) {

    size_t nPoints = 1;
    size_t dim     = 1;

    auto result = MonteCarlo::TransposerIndices ( nPoints, dim );

    ASSERT_EQ ( result.size(), nPoints * dim );

    EXPECT_EQ ( result[0], 0 );

}

TEST ( TransposerIndices, MultipleIntervalSingleSample ) {

    size_t nPoints = 4;
    size_t dim     = 1;

    auto result = MonteCarlo::TransposerIndices ( nPoints, dim );

    ASSERT_EQ ( result.size(), nPoints * dim );

    std::vector<size_t> expected {

        0, 1, 2, 3 

    };

    for ( auto i = 0; i < result.size(); i++ ) {

        EXPECT_EQ ( result[i] , expected[i] );

    }

}

TEST ( TransposerIndices, SingleIntervalMultipleSamples ) {

    size_t nPoints = 1;
    size_t dim     = 5;

    auto result = MonteCarlo::TransposerIndices ( nPoints, dim );

    ASSERT_EQ ( result.size(), nPoints * dim );

    std::vector<size_t> expected {

        0, 1, 2, 3, 4 

    };

    for ( auto i = 0; i < result.size(); i++ ) {

        EXPECT_EQ ( result[i] , expected[i] );

    }

}

TEST ( TransposerIndices, MultipleIntervalMultipleSamples ) {

    size_t nPoints = 3;
    size_t dim     = 5;

    auto result = MonteCarlo::TransposerIndices ( nPoints, dim );

    ASSERT_EQ ( result.size(), nPoints * dim );

    std::vector<size_t> expected {

         0,  3,  6,  9, 12, 
         1,  4,  7, 10, 13, 
         2,  5,  8, 11, 14 

    };

    for ( auto i = 0; i < result.size(); i++ ) {

        EXPECT_EQ ( result[i] , expected[i] );

    }

}

TEST ( TransposerIndices, WrongNPoints ) {

    size_t nPoints = 0;
    size_t dim     = 2;

    EXPECT_THROW ({

        try {

            auto result = MonteCarlo::TransposerIndices ( nPoints, dim );

        } catch ( const std::exception& e ) {

            EXPECT_STREQ (
                "TransposerIndices: num of points must be greater than unity", 
                e.what()
            );

            throw;

        }


    }, std::runtime_error );

}

TEST ( TransposerIndices, WrongDim ) {

    size_t nPoints = 5;
    size_t dim     = 0;

    EXPECT_THROW ({

        try {

            auto result = MonteCarlo::TransposerIndices ( nPoints, dim );

        } catch ( const std::exception& e ) {

            EXPECT_STREQ (
                "TransposerIndices: dimension must be greater than unity", 
                e.what()
            );

            throw;

        }


    }, std::runtime_error );

}

TEST ( LatinHypercubeSampling, SingleIntervalSingleSample ) {

    size_t nPoints = 1;
    size_t dim     = 1;

    auto result = 
        MonteCarlo::LHS < size_t, double > ( nPoints, dim );

    double range = 1.0 / nPoints;

    ASSERT_EQ ( result.size(), nPoints * dim );

    for ( auto i = 0; i < dim; i++ ) {

        std::vector <bool> indices ( nPoints, false );

        for ( auto j = 0; j < nPoints; j++ ) {
            
            for ( auto k = 0; k < nPoints; k++ ) {

                if (
                    result[i + j * dim] >=  k    * range &&
                    result[i + j * dim] <= (k+1) * range
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

TEST ( LatinHypercubeSampling, MultipleIntervalsSingleSample ) {

    size_t nPoints = 40;
    size_t dim     = 1;

    auto result = 
        MonteCarlo::LHS < size_t, double > ( nPoints, dim );

    double range = 1.0 / nPoints;

    ASSERT_EQ ( result.size(), nPoints * dim );

    for ( auto i = 0; i < dim; i++ ) {

        std::vector <bool> indices ( nPoints, false );

        for ( auto j = 0; j < nPoints; j++ ) {
            
            for ( auto k = 0; k < nPoints; k++ ) {

                if (
                    result[i + j * dim] >=  k    * range &&
                    result[i + j * dim] <= (k+1) * range
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

    size_t nPoints = 1;
    size_t dim     = 50;

    auto result = 
        MonteCarlo::LHS < size_t, double > ( nPoints, dim );

    double range = 1.0 / nPoints;

    ASSERT_EQ ( result.size(), nPoints * dim );

    for ( auto i = 0; i < dim; i++ ) {

        std::vector <bool> indices ( nPoints, false );

        for ( auto j = 0; j < nPoints; j++ ) {
            
            for ( auto k = 0; k < nPoints; k++ ) {

                if (
                    result[i + j * dim] >=  k    * range &&
                    result[i + j * dim] <= (k+1) * range
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


TEST ( LatinHypercubeSampling, MultipleIntervalsMultipleSamples ) {

    size_t nPoints = 40;
    size_t dim     = 50;

    auto result = 
        MonteCarlo::LHS < size_t, double > ( nPoints, dim );

    double range = 1.0 / nPoints;

    ASSERT_EQ ( result.size(), nPoints * dim );

    for ( auto i = 0; i < dim; i++ ) {

        std::vector <bool> indices ( nPoints, false );

        for ( auto j = 0; j < nPoints; j++ ) {
            
            for ( auto k = 0; k < nPoints; k++ ) {

                if (
                    result[i + j * dim] >=  k    * range &&
                    result[i + j * dim] <= (k+1) * range
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


