/**
  * @file main.cpp
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "SurrogateModel.hpp" 
#include <iostream> 

int main () {

    typedef double Float; 
    typedef std::complex<Float> Complex; 

    typedef std::vector<Float>   VectorF;
    typedef std::vector<Complex> VectorC;

    typedef Analytical::MassSpringDamper<size_t,Float,Complex> AnalyticalModel;
    typedef Surrogate::IntrusivePCE SurrogateModel;


    // Define analytical model 

    VectorF Masses { 2.0, 2.0 };
    VectorF Damper { 5.0, 2.5 };
    VectorF Spring { 20.0, 10.0 };

    AnalyticalModel SDModel ( Masses, Damper, Spring );


    // Define deterministic load 

    VectorC Force { 1.0, 1.0 };
    Float   Omega = 0.0; 


    // Create and train surrogate model 

    SurrogateModel NPCE ( &SDModel, Omega );

    NPCE.SetIndices (
        10,  // max index 
        10   // max sum of indices 
    );

    NPCE.Train ( Force );


    // Check with two points 

    VectorC RandomParts { 

        0.0, 0.0  // point #1
        // 0.0, 2.0   // point #2

    };

    NPCE.PrintCoeffs();

    auto dispResult = NPCE.ComputeResponse ( RandomParts );

    for ( auto i = 0; i < dispResult.size(); i++ ) {
        std::cout << dispResult[i] << std::endl;
    }

} // main 

