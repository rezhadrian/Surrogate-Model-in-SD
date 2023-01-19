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
  * @file Truncations.hpp
  *
  * @brief declare functions to truncate indices of orthonormal polynomials.
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef TRUNCATIONS_DECLARATIONS 
#define TRUNCATIONS_DECLARATIONS 

#include "LibrariesLoader.hpp" 

namespace BasisFunctions {

    // Some aliases
    template < typename T >
    using Vector = std::vector <T>;


    template < typename T >
    /**
      * Remove tupples with sum larger than SumPMax.
      * 
      * @tparam T a type of integer e.g. size_t.
      *
      * @param MultiIndex tupples of indices of orthonormal polynomials.
      * @param        dim number of elements inside each tupple.
      * @param    SumPMax largest allowable sum of indices in each tupple.
      */
    void TotalTruncation ( Vector<T>& MultiIndex, T dim, T SumPMax );


} // BasisFunctions 

#ifndef TRUNCATIONS_IMPLEMENTATIONS 
    #include "Truncations_imp.hpp" 
#endif 

#endif // TRUNCATIONS_DECLARATIONS 

