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

    typedef Surrogate::IntrusivePCE  IntrusivePCE;
    typedef Surrogate::IntrusiveRPCE IntrusiveRPCE;

    typedef Surrogate::NonIntrusivePCE NonIntrusivePCE;
    typedef Surrogate::NonIntrusiveRPCE NonIntrusiveRPCE;


    // =======================================================================
    // Create analytical model and define load 
    // =======================================================================

    VectorF Masses { 1.0, 1.0 };
    VectorF Damper { 1.0, 0.5 };
    VectorF Spring { 20.0, 10.0 };

    AnalyticalModel SDModel ( Masses, Damper, Spring );

    VectorC Force { 1.0, 1.0 };
    Float   Omega = 2.2; 


    // =======================================================================
    // Create and train intrusive PCE surrogate model 
    // =======================================================================

    IntrusivePCE iPCE ( &SDModel, Omega );

    iPCE.SetIndices (
        8,  // max index 
        8   // max sum of indices 
    );

    iPCE.Train ( Force );


    // =======================================================================
    // Create and train intrusive RPCE surrogate model 
    // =======================================================================

    IntrusiveRPCE iRPCE ( &SDModel, Omega );

    iRPCE.SetNumIndices ( 1, 4 );
    iRPCE.SetDenIndices ( 3, 4 );

    iRPCE.Train ( Force );


    // =======================================================================
    // Create and train non-intrusive PCE surrogate model 
    // =======================================================================

    NonIntrusivePCE niPCE ( &SDModel, Omega );

    niPCE.SetIndices ( 8, 8 ); 

    auto TrainSet = MonteCarlo::LHS<size_t,Float>( 80, 2 );

    MonteCarlo::ConvertLHStoStdNorm<size_t,Float>(TrainSet);

    niPCE.Train ( Force, TrainSet );


    // =======================================================================
    // Create and train non-intrusive RPCE surrogate model 
    // =======================================================================

    NonIntrusiveRPCE niRPCE ( &SDModel, Omega );

    niRPCE.SetNumIndices ( 1, 4 );
    niRPCE.SetDenIndices ( 3, 4 );

    niRPCE.Train ( Force, TrainSet );


    // =======================================================================
    // Generate random input and store in csv format 
    // =======================================================================

    VectorC RandomInput = 
        MonteCarlo::RandomSampling<size_t,Float,Complex> ( 10000, 2 );

    VectorC trainInput ( TrainSet.size() );

    std::transform ( 

        TrainSet.begin(), 
        TrainSet.end(), 

        trainInput.begin(), 

        [](const auto m) { return Complex(m); }

    );

    std::ofstream randomInput ( "input.csv" );

    for ( auto i = 0; i < RandomInput.size(); i++ ) {

        randomInput << RandomInput[i] << "\n";

    }

    randomInput.close();


    // =======================================================================
    // Perform Monte Carlo simulation on intrusive PCE model 
    // =======================================================================

    auto iPCEoutput = iPCE.ComputeResponse ( RandomInput );

    std::ofstream iPCEstream ( "iPCE.csv" );

    for ( auto i = 0; i < iPCEoutput.size(); i++ ) {

        iPCEstream << iPCEoutput[i] << "\n";

    }

    iPCEstream.close();


    // =======================================================================
    // Perform Monte Carlo simulation on intrusive RPCE model 
    // =======================================================================

    auto iRPCEoutput = iRPCE.ComputeResponse ( RandomInput );

    std::ofstream iRPCEstream ( "iRPCE.csv" );

    for ( auto i = 0; i < iRPCEoutput.size(); i++ ) {

        iRPCEstream << iRPCEoutput[i] << "\n";

    }

    iRPCEstream.close();


    // =======================================================================
    // Perform Monte Carlo simulation on non-intrusive PCE model 
    // =======================================================================

    auto niPCEoutput = niPCE.ComputeResponse ( RandomInput );

    std::ofstream niPCEstream ( "niPCE.csv" );

    for ( auto i = 0; i < niPCEoutput.size(); i++ ) {

        niPCEstream << niPCEoutput[i] << "\n";

    }

    niPCEstream.close();


    // =======================================================================
    // Perform Monte Carlo simulation on non-intrusive PCE model 
    // =======================================================================


    auto niRPCEoutput = niRPCE.ComputeResponse ( RandomInput );

    std::ofstream niRPCEstream ( "niRPCE.csv" );

    for ( auto i = 0; i < niRPCEoutput.size(); i++ ) {

        niRPCEstream << niRPCEoutput[i] << "\n";

    }

    niRPCEstream.close();

} // main 

