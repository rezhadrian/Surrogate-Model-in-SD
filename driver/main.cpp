/**
  * @file main.cpp
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "SurrogateModel.hpp" 
#include <iostream> 
#include <fstream> 

int main () {

    typedef double Float; 
    typedef std::complex<Float> Complex; 

    typedef std::vector<Float>   VectorF;
    typedef std::vector<Complex> VectorC;

    typedef Analytical::MassSpringDamper<size_t,Float,Complex> AnalyticalModel;
    typedef Surrogate::IntrusivePCE SurrogateModel;
    typedef Surrogate::IntrusiveRPCE SurrogateModelR;


    // Define analytical model 

    VectorF Masses { 1.0, 1.0 };
    VectorF Damper { 1.0, 0.5 };
    VectorF Spring { 20.0, 10.0 };

    AnalyticalModel SDModel ( Masses, Damper, Spring );


    // Define deterministic load 

    VectorC Force { 1.0, 1.0 };
    Float   Omega = 2.2; 


    // Create and train surrogate model 

    SurrogateModel NPCE ( &SDModel, Omega );

    NPCE.SetIndices (
        4,  // max index 
        4   // max sum of indices 
    );

    NPCE.Train ( Force );

    SurrogateModelR NRPCE ( &SDModel, Omega );

    NRPCE.SetNumIndices ( 1, 4 );
    NRPCE.SetDenIndices ( 2, 4 );

    NRPCE.Train ( Force );


    // Check with two points 

    VectorC RandomParts = 
        MonteCarlo::RandomSampling<size_t,Float,Complex> ( 10000, 2 );

    std::ofstream input ( "input.csv" );

    for ( auto i = 0; i < RandomParts.size(); i++ ) {

        input << RandomParts[i] << "\n";

    }

    input.close();

    // VectorC RandomParts { 
    //
    //     0.0, 0.0  // point #1
    //     // 0.0, 2.0   // point #2
    //
    // };
    //
    // NPCE.PrintCoeffs();

    auto dispResult = NPCE.ComputeResponse ( RandomParts );

    auto anotherResult = NRPCE.ComputeResponse ( RandomParts );

    // for ( auto i = 0; i < anotherResult.size(); i++ ) {
    //
    //     std::cout << anotherResult[i] << std::endl;
    //
    // }

    std::ofstream output ( "output.csv" );

    for ( auto i = 0; i < dispResult.size(); i++ ) {

        output << dispResult[i] << "\n";

    }

    output.close();

    std::ofstream output2 ( "output2.csv" );

    for ( auto i = 0; i < anotherResult.size(); i++ ) {

        output2 << anotherResult[i] << "\n";

    }

    output2.close();
    // for ( auto i = 0; i < dispResult.size(); i++ ) {
    //     std::cout << dispResult[i] << std::endl;
    // }

} // main 

