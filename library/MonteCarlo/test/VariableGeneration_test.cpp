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

TEST ( ConvertLHStoStdNorm, MockLHSResult ) {

    std::vector<double> LHSResult {

        0.5000, 0.5250, 0.5500, 0.5750, 0.6000, 0.6250, 0.6500, 0.6750, 
        0.7000, 0.7250, 0.7500, 0.7750, 0.8000, 0.8250, 0.8500, 0.8750, 
        0.5000, 0.4750, 0.4500, 0.4250, 0.4000, 0.3750, 0.3500, 0.3250,
        0.3000, 0.2750, 0.2500, 0.2250, 0.2000, 0.1750, 0.1500, 0.1250,
        0.9210, 0.9215, 0.9250, 0.9300, 0.9350, 0.9400, 0.9450, 0.9500,
        0.9550, 0.9600, 0.9650, 0.9700, 0.9750, 0.9800, 0.9850, 0.9900,
        0.9950, 0.9960, 0.9970, 0.9980, 0.9990, 0.9995, 0.9996, 0.9997,
        0.9998, 0.9999, 0.0790, 0.0785, 0.0750, 0.0700, 0.0650, 0.0600,
        0.0550, 0.0500, 0.0450, 0.0400, 0.0350, 0.0300, 0.0250, 0.0200,
        0.0150, 0.0100, 0.0050, 0.0040, 0.0030, 0.0020, 0.0010, 0.0005,
        0.0004, 0.0003, 0.0002, 0.0001

    };

    std::vector<double> StdNorm {

         0.000000,  0.062707,  0.125661,  0.189118,  0.253347,  0.318639,
         0.385320,  0.453762,  0.524401,  0.597760,  0.674490,  0.755415,
         0.841621,  0.934589,  1.036433,  1.150349,  0.000000, -0.062707,
        -0.125661, -0.189118, -0.253347, -0.318639, -0.385320, -0.453762,
        -0.524401, -0.597760, -0.674490, -0.755415, -0.841621, -0.934589,
        -1.036433, -1.150349,  1.411830,  1.415234,  1.439531,  1.475791,
         1.514102,  1.554774,  1.598193,  1.644854,  1.695398,  1.750686,
         1.811911,  1.880794,  1.959964,  2.053749,  2.170090,  2.326348,
         2.575829,  2.652070,  2.747781,  2.878162,  3.090232,  3.290527,
         3.352795,  3.431614,  3.540084,  3.719016, -1.411830, -1.415234,
        -1.439531, -1.475791, -1.514102, -1.554774, -1.598193, -1.644854,
        -1.695398, -1.750686, -1.811911, -1.880794, -1.959964, -2.053749,
        -2.170090, -2.326348, -2.575829, -2.652070, -2.747781, -2.878162,
        -3.090232, -3.290527, -3.352795, -3.431614, -3.540084, -3.719016

    };

    MonteCarlo::ConvertLHStoStdNorm<size_t,double> ( LHSResult );

    double tol = 0.000001;

    for ( auto i = 0; i < LHSResult.size(); i++ ) {

        EXPECT_NEAR ( LHSResult[i], StdNorm[i], tol);
        
    }
}

// TEST ( GenerateRVs, MockStdNorm ) {
//
//     typedef MonteCarlo::DenseSymMatrix <size_t, double> SM;
//     typedef MonteCarlo::DenseLTriangularMatrix <size_t,double> LT;
//     typedef std::vector<double> Vector;
//
//     size_t dim = 4;
//
//     Vector StdNormRVs {
//
//          0.1,  0.2,  0.3,  0.4, 
//          0.4,  0.3,  0.2,  0.1, 
//          0.5,  0.5,  0.4,  0.4, 
//          0.3,  0.3,  0.2,  0.2, 
//         -0.1, -0.2, -0.3, -0.4, 
//         -0.4, -0.3, -0.2, -0.1, 
//         -0.5, -0.5, -0.4, -0.4, 
//         -0.3, -0.3, -0.2, -0.2 
//
//     };
//
//
//     SM Corr ( dim );
//
//     Corr ( 0, 0 ) =  1.0000; Corr ( 1, 0 ) = -0.7601; 
//     Corr ( 2, 0 ) =  0.9343; Corr ( 3, 0 ) =  0.9691; 
//     Corr ( 1, 1 ) =  1.0000; Corr ( 2, 1 ) = -0.6293; 
//     Corr ( 3, 1 ) = -0.7881; Corr ( 2, 2 ) =  1.0000; 
//     Corr ( 3, 2 ) =  0.8920; Corr ( 3, 3 ) =  1.0000; 
//
//     std::vector<std::function<double(double)>> ICDFs;
//
//
//     ICDFs.emplace_back (
//
//         [](const auto m) {
//             return MonteCarlo::InvLogNormCDF<double> ( m, 0.0, 1.0 );
//         }
//
//     );
//
//     ICDFs.emplace_back (
//
//         [](const auto m) {
//             return MonteCarlo::InvLogNormCDF<double> ( m, 0.1, 0.5 );
//         }
//
//     );
//
//     ICDFs.emplace_back (
//
//         [](const auto m) {
//             return MonteCarlo::InvLogNormCDF<double> ( m, 0.3, 0.2 );
//         }
//
//     );
//
//     ICDFs.emplace_back (
//
//         [](const auto m) {
//             return MonteCarlo::InvLogNormCDF<double> ( m, 0.3, 0.5 );
//         }
//
//     );
//
//     Vector result = 
//         MonteCarlo::GenerateRVs < double, LT, std::function<double(double)>
//
//         > ( StdNormRVs, Corr, ICDFs, dim );
//
//     Vector expected {
//
//         1.105170918075648, 1.135389327255325, 1.410165929127596, 
//         1.470494476447416, 1.491824697641270, 1.046499638624557, 
//         1.485234831925537, 1.636457192059457, 1.648721270700128, 
//         1.075113813834130, 1.541258060148288, 1.762968747246328, 
//         1.349858807576003, 1.087037295008275, 1.457739423492359, 
//         1.577350001291581, 0.904837418035960, 1.075756772448067, 
//         1.292130778906117, 1.239119785606122, 0.670320046035639, 
//         1.167131562286532, 1.226822022500117, 1.113453385295951, 
//         0.606530659712634, 1.136068332899878, 1.182228237765192, 
//         1.033551390650895, 0.740818220681718, 1.123607040686559, 
//         1.249961941775021, 1.155177226930297
//         
//     };
//
//     double tol = 0.000000000000001;
//
//     ASSERT_EQ ( result.size(), expected.size() );
//
//     for ( auto i = 0; i < result.size(); i++ ) {
//
//         EXPECT_NEAR ( result[i], expected[i], tol );
//
//     }
//
// }
//
