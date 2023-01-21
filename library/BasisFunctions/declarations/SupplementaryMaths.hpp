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
  * @file SupplementaryMaths_imp.hpp
  *
  * @brief implement additional maths feature not implemented in <cmath> 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef SUPPLEMENTARY_MATHS_DECLARATIONS 
#define SUPPLEMENTARY_MATHS_DECLARATIONS 

#include "LibrariesLoader_BF.hpp" 

namespace BasisFunctions {


    template < typename Z > 
    /**
      * Compute n! 
      * @tparam Z a type of non-negative integer e.g. size_t 
      */
    Z Factorial ( const Z n ); 


    template < typename Z, typename C >
    /**
      * Overload multiplication for integer and complex operand 
      *
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam C a class comparable with std::complex 
      */
    C operator* ( const Z integer, const C complex );


    template < typename Z >
    /**
      * Compute binomial coefficient nCk recursively.
      *
      * @tparam Z a type of non-negative integer e.g. size_t, int.
      */
    Z Binomial ( const Z n, const Z k );


} // BasisFunctions 

#ifndef SUPPLEMENTARY_MATHS_IMPLEMENTATIONS 
    #include "SupplementaryMaths_imp.hpp" 
#endif 

#endif // SUPPLEMENTARY_MATHS_DECLARATIONS 

