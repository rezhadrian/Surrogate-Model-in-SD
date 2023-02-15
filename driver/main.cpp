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
    VectorF Damper { 1.0, 0.5 };
    VectorF Spring { 2.0, 1.0 };

    AnalyticalModel SDModel ( Masses, Damper, Spring );


    // Define deterministic load 

    VectorC Force { 1.0, 1.0 };
    Float   Omega = 0.0; 


    // Create and train surrogate model 

    SurrogateModel NPCE ( &SDModel, Omega );

    NPCE.SetIndices (
        10, 10
    );

    NPCE.Train ( Force );


    // Check with two points 

    VectorC RandomParts { 

        1.0, 0.0  // point #1
        // 0.0, 2.0   // point #2

    };

    auto dispResult = NPCE.ComputeResponse ( RandomParts );

    for ( auto i = 0; i < dispResult.size(); i++ ) {
        std::cout << dispResult[i] << std::endl;
    }

} // main 

