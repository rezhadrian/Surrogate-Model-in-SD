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


// Implemented in SupplementaryMaths_imp_BF.hpp 

namespace BasisFunctions {


    template < typename Z > 
    /**
      * Compute n! recursively 
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
      * Compute binomial coefficient nCk recursively 
      *
      * @tparam Z a type of non-negative integer e.g. size_t 
      */
    Z Binomial ( const Z n, const Z k );


} // BasisFunctions : SupplementaryMaths 


// Implemented in MultiIndex_imp.hpp 

namespace BasisFunctions {

    using Marker = std::vector <bool>;


    template < typename Z >
    /**
      * Generate Indices of orthonormal polynomials 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      *
      * @param   dim number of variables 
      * @param  pMax maximum order of the polynomial 
      * @return sequence of polynomial indices 
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
    void MultiIndexRecursive ( 
        const Z nSet, Vector<Z>& Subset, Marker& Active 
    );

    
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


    template < typename Z, typename C >
    /**
      * Evaluate tensor products of NORMALIZED Hermite polynomials.
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


    template < typename Z, typename C >
    /**
      * Evaluate probabilist Hermite polynomial of given index at given position.
      * The polynomial is NOT normalized.
      *
      * @tparam Z a type of non-negative integer e.g. size_t
      * @tparam C a type of floating complex number e.g. std::complex<float> 
      * 
      * @param  index indicates which polynomial to use 
      * @param      x argument to evaluate the polynomial 
      * @return H_{idx} ( x ) 
      */
    C HermitePolynomial ( const Z index, const C x );


} // BasisFunctions : HermitePolynomials 


#ifndef SUPPLEMENTARY_MATHS_IMPLEMENTATIONS 
    #include "SupplementaryMaths_imp_BF.hpp" 
#endif 

#ifndef MULTI_INDEX_IMPLEMENTATIONS 
    #include "MultiIndex_imp.hpp" 
#endif 

#ifndef TRUNCATIONS_IMPLEMENTATIONS 
    #include "Truncations_imp.hpp" 
#endif 

#ifndef HERMITE_POLYNOMIALS_IMPLEMENTATIONS 
    #include "HermitePolynomials_imp.hpp" 
#endif 

#endif // BASIS_FUNCTIONS_DECLARATIONS 

