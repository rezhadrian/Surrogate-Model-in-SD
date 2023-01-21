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
  * @file MultiIndex.hpp
  *
  * @brief declare functions to calculate indices of orthonormal polynomials.
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef MULTI_INDEX_DECLARATIONS 
#define MULTI_INDEX_DECLARATIONS 

#include "LibrariesLoader_BF.hpp" 
#include "SupplementaryMaths.hpp" 

namespace BasisFunctions {

    // Some Aliases
    template < typename T >
    using Vector = std::vector <T>;
    using Marker = std::vector <bool>;


    template < typename Z >
    /**
      * Generate Indices of orthonormal polynomials.
      * 
      * @tparam Z a type of non-negative integer e.g. size_t.
      *
      * @param   dim number of variables.
      * @param  pMax maximum order of the polynomial.
      * @return sequence of polynomial indices.
      */
    Vector<Z> MultiIndex ( const Z dim, const Z pMax );


    template < typename Z >
    /**
      * Generate a unique tupple of indices for MultiIndex function. 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t.
      *
      * @param   nSet the tupple is generated from the first nSet integers.
      * @param Subset save the tupple here. Also contains all previous tupples.
      * @param Active track all previous tupples that have been saved.
      */
    void MultiIndexRecursive ( const Z nSet, Vector<Z>& Subset, Marker& Active );


} // BasisFunctions 

#ifndef MULTI_INDEX_IMPLEMENTATIONS 
    #include "MultiIndex_imp.hpp" 
#endif 

#endif // MULTI_INDEX_DECLARATIONS 

