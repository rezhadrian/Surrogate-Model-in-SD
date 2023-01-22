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
  * @file InverseCDF.hpp
  *
  * @brief declare functions to compute inverse CDF of standard normal dist
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef INVERSE_CDF_DECLARATIONS 
#define INVERSE_CDF_DECLARATIONS 

#include "LibrariesLoader_MC.hpp" 
#include "SupplementaryMaths_MC.hpp" 

namespace MonteCarlo {

    template < typename T >
    using Vector = std::vector <T>;


    template < typename Z, typename R >
    /**
      * Apply Inverse CDF function to a vector of uniformly distributed RVs
      * Results overwrite inputs 
      *
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      *
      * @param UniformSamples vector with elements in (0,1) 
      */
    void InverseCDFs ( Vector<R>& UniformSamples );


    template < typename Z, typename R >
    /**
      * Compute inverse standard normal CDF of a given quantile 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      * 
      * @return z such that CDF(z) = quantile 
      */
    R InverseSingleCDF ( const R quantile );


} // MonteCarlo 

#ifndef INVERSE_CDF_IMPLEMENTATIONS 
    #include "InverseCDF_imp.hpp" 
#endif 

#endif // INVERSE_CDF_DECLARATIONS 

