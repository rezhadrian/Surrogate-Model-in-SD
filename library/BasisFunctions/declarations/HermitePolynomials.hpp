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
  * @file HermitePolynomials.hpp
  *
  * @brief declare functions to calculate probabilist Hermite polynomials.
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef HERMITE_POLYNOMIALS_DECLARATIONS 
#define HERMITE_POLYNOMIALS_DECLARATIONS 

#include "LibrariesLoader.hpp" 

namespace BasisFunctions {

    // Some Aliases 
    template < typename T >
    using Vector = std::vector<T>;


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
    Vector<C> HermitePolynomials ( Vector<Z>& indices, Vector<C>& X, Z dim );


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
    C HermitePolynomial ( Z index, C x );

    
    template < typename Z, typename C >
    /**
      * Compute power function recursively instead of using std::pow.
      * 
      * @tparam Z a type of non-negative integer e.g. size_t
      * @tparam C a type of floating complex number e.g. std::complex<float> 
      * @return x^{power} 
      */
    C Pow ( Z power, C x );

    template < typename Z >
    /**
      * Compute factorial function recursively.
      *
      * @tparam Z a type of non-negative integer e.g. size_t
      * @return n!
      */
    Z Factorial ( Z n );

    template < typename Z, typename C >
    C Multiply ( Z index, C x );

} // BasisFunctions 

#ifndef HERMITE_POLYNOMIALS_IMPLEMENTATIONS 
    #include "HermitePolynomials_imp.hpp" 
#endif 

#endif // HERMITE_POLYNOMIALS_DECLARATIONS 

