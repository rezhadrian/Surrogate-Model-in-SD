/**
  * @file LatinHypercubeSampling_imp.hpp
  *
  * @brief implement functions to produce latin hypercube sample
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef LATIN_HYPERCUBE_SAMPLING_IMPLEMENTATIONS 
#define LATIN_HYPERCUBE_SAMPLING_IMPLEMENTATIONS 

#ifndef MONTE_CARLO_DECLARATIONS 
    #include "MonteCarlo.hpp" 
#endif 


namespace MonteCarlo {

    template < 

        typename Z, typename R,
               
        class RandomDevice,
        class RandomEngine,
        class Distribution,
        class Shuffler 

    >
    Vector<R> LHS ( const Z nInterval, const Z nSample ) {

        if ( nInterval < 1 ) {

            throw std::runtime_error (
                "LHS: interval must be positive integer"
            );

        }

        if ( nSample < 1 ) {

            throw std::runtime_error (
                "LHS: number of sample must be positive integer"
            );

        }

        R range = 1.0 / nInterval;

        RandomDevice device;

        RandomEngine generator ( device () );
        Shuffler     shuffler  ( device () );

        Distribution RandomVariable ( 0.0, range );


        Vector<Z> IntervalIndices ( nInterval );

        std::iota (

            IntervalIndices.begin(),
            IntervalIndices.end(),
            0

        );

        Vector<R> result ( nInterval * nSample );

        for ( auto i = 0; i < nSample; i++ ) {

            std::shuffle (

                IntervalIndices.begin(),
                IntervalIndices.end(),
                shuffler

            );

            std::transform (

                IntervalIndices.begin(),
                IntervalIndices.end(),

                result.begin() + i * nInterval,

                [range, &RandomVariable, &generator ]( auto m ) {
                    
                    return m * range + RandomVariable ( generator );
                    
                }

            );

        }

        return result;

    }

} // MonteCarlo : LHS 


#endif // LATIN_HYPERCUBE_SAMPLING_IMPLEMENTATIONS 

