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
  * @file BasisFunctions.hpp 
  *
  * @brief declarations of all functions needed to generate basis functions 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef BASIS_FUNCTIONS_DECLARATIONS 
#define BASIS_FUNCTIONS_DECLARATIONS 

#include "LibrariesLoader_BF.hpp" 


template < typename T >
using Vector = std::vector <T>;


// Implemented in MultiIndex_imp.hpp 

namespace BasisFunctions {

    using Marker = std::vector <bool>;


    template < typename Z >
    /**
      * Generate Indices of orthonormal polynomials 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      *
      * @param  dim number of variables 
      * @param  pMax maximum order of the polynomial 
      * @return sequence of polynomial indices 
      */
    Vector<Z> MultiIndex ( const Z dim, const Z pMax );

    
} // BasisFunctions : MultiIndex 


// Implemented in Truncations_imp.hpp 

namespace BasisFunctions {


    template < typename Z >
    /**
      * Remove tupples with sum larger than SumPMax 
      * 
      * @tparam T a type of integer e.g. size_t 
      *
      * @param MultiIndex tupples of indices of orthonormal polynomials 
      * @param        dim number of elements inside each tupple 
      * @param    SumPMax largest allowable sum of indices in each tupple 
      */
    void TotalTruncation ( 
        Vector<Z>& MultiIndex, const Z dim, const Z SumPMax 
    );


} // BasisFunctions : Truncations 


// Implemented in HermitePolynomials_imp.hpp 

namespace BasisFunctions {


    template < typename Z, typename R, typename C >
    /**
      * Evaluate tensor prod. of NORMALIZED probabilist Hermite polynomials.
      *
      * @tparam Z a type of non-negative integer e.g. size_t
      * @tparam C a type of floating complex number e.g. std::complex<float> 
      *
      * @param  indices tupples of indices for the polynomial 
      * @param        X vector of arguments to evaluate the polynomial
      * @param      dim define how many polynomials are multiplied together 
      * @return vector of tensor products of Hermite polynomials
      */
    Vector<C> HermitePolynomials (
        const Vector<Z>& indices, const Vector<C>& X, const Z dim 
    );


} // BasisFunctions : HermitePolynomials 


// Implemented in TripleHermite_imp.hpp 

namespace BasisFunctions {

    template < typename Z, typename R > 
    /**
      * Compute expected value of products of three Hermite polynomials 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating number e.g. double 
      * 
      * @return E( H_i, H_j, H_k ) 
      */
    R ETripleHermite ( const Z i, const Z j, const Z k );


} // BasisFunctions : TripleHermite 


#ifndef MULTI_INDEX_IMPLEMENTATIONS 
    #include "MultiIndex_imp.hpp" 
#endif 

#ifndef TRUNCATIONS_IMPLEMENTATIONS 
    #include "Truncations_imp.hpp" 
#endif 

#ifndef HERMITE_POLYNOMIALS_IMPLEMENTATIONS 
    #include "HermitePolynomials_imp.hpp" 
#endif 

#ifndef TRIPLE_HERMITE_IMPLEMENTATIONS 
    #include "TripleHermite_imp.hpp" 
#endif 

#endif // BASIS_FUNCTIONS_DECLARATIONS 

