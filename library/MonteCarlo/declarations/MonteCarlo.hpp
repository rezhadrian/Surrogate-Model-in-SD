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
  * @file MonteCarlo.hpp
  *
  * @brief declare functions required to run computer experiments 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef MONTE_CARLO_DECLARATIONS 
#define MONTE_CARLO_DECLARATIONS 

#include "LibrariesLoader_MC.hpp" 


#ifndef STATISTIC_DIST_DECLARATIONS 
    #include "StatisticDist.hpp" 
#endif 

#ifndef LINALG_DECLARATIONS 
    #include "LinAlg.hpp" 
#endif 


// Implemented in LatinHypercubeSampling_imp.hpp 

namespace MonteCarlo {


    template < typename Z, typename R,
               
        class RandomDevice,   // e.g. std::random_device 
        class RandomEngine,   // e.g. std::default_random_engine 
        class Distribution,   // e.g. std::uniform_real_distribution 
        class Shuffler        // e.g. std::mt19937 

    >
    /**
      * Perform latin hypercube sampling with given sizes 
      *
      * @tparam Z a type of non=negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      * 
      * @param nInterval number of intervals in [0,1] 
      * @param nSample   number of sampled [0,1] 
      * @return vector { Sample1, Sample2, ... } 
      */
    Vector<R> LHS ( const Z nInterval, const Z nSample );


} // MonteCarlo : LatinHypercubeSampling 


// Implemented in VariableGeneration_imp.hpp 

namespace MonteCarlo {


    template < typename Z, typename R >
    /**
      * Convert LHS result to sets of standard normal random variables 
      *
      * @tparam Z a type of non=negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      */
    void ConvertLHS ( Vector<R>& LHSResult );


    template < typename Z, typename R, 

        class SymMatrix,         // e.g. MonteCarlo::DenseSymMatrix 
        class LTriangularMatrix, // e.g. MonteCarlo::DenseLTriangularMatrix 
        class Function 

    >
    /**
      * Convert standard normal RVs to correlated RVs of other distribution 
      * 
      * @tparam Z a type of non=negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      */
    Vector<R> ConvertStdNorm ( 

        const Vector<R>& StdNormRVs, 
        const SymMatrix& Correlation, 
        const Function & ICDF,
        const Z dim

    );


} // MonteCarlo : ConvertLHS 


#ifndef LATIN_HYPERCUBE_SAMPLING_IMPLEMENTATIONS 
    #include "LatinHypercubeSampling_imp.hpp" 
#endif 

#ifndef VARIABLE_GENERATION_IMPLEMENTATIONS 
    #include "VariableGeneration_imp.hpp" 
#endif 

#endif // MONTE_CARLO_DECLARATIONS 

