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
  * @file LatinHypercubeSampling.hpp
  *
  * @brief declare functions to produce latin hypercube sample
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef LATIN_HYPERCUBE_SAMPLING_DECLARATIONS 
#define LATIN_HYPERCUBE_SAMPLING_DECLARATIONS 

#include "LibrariesLoader_MC.hpp" 

namespace MonteCarlo {

    template < typename T >
    using Vector = std::vector <T>;


    template < typename Z, typename R,
               
               class RandomDevice,
               class RandomEngine,
               class Distribution,
               class Shuffler 

    >
    /**
      * Perform latin hypercube sampling with given sizes 
      *
      * @tparam Z a type of non=negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      * 
      * @tparam RandomDevice non-deterministic random number generator
      *         e.g. std::random_device 
      * @tparam RandomEngine pseudo-random number generator between 0 and 1
      *         e.g. std::default_random_engine 
      * @tparam Distribution distribution within each interval 
      *         e.g. std::uniform_real_distribution<double> 
      * @tparam Shuffler pseudo random shuffler 
      *         e.g. std::mt19937
      * 
      * @param nInterval number of intervals in [0,1] 
      * @param nSample   number of sampled [0,1] 
      * @return vector { Sample1, Sample2, ... } 
      */
    Vector<R> LHS ( const Z nInterval, const Z nSample );


} // MonteCarlo 

#ifndef LATIN_HYPERCUBE_SAMPLING_IMPLEMENTATIONS 
    #include "LatinHypercubeSampling_imp.hpp" 
#endif 

#endif // LATIN_HYPERCUBE_SAMPLING_DECLARATIONS 

