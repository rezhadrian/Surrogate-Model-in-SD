/**
  * This file is part of a Surrogate Model library for structural dynamics. 
  *
  * The library is free: you can distribute it and/or modify it.
  * The library is distributed in the hope that it can be useful and helpful,
  * particularly for students and self - learning individuals.
  * 
  * The library comes WITHOUT ANY WARRANTY: without even the implied warranty 
  * of MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE.
  *
  * The library is licensed under the GNU General Public License v3.0
  *
  *
  * The library is based on works by:
  * - Felix Schneider
  * - Iason Papaioannou
  * 
  * Chair of Structural Mechanics 
  * Engineering Risk Analysis Group 
  *
  * Technische Universitaet Muenchen
  *
  * www.cee.ed.tum.de/bm
  * www.cee.ed.tum.de/era 
  */ 

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

#ifndef LATIN_HYPERCUBE_SAMPLING_DECLARATIONS 
    #include "LatinHypercubeSampling.hpp" 
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

