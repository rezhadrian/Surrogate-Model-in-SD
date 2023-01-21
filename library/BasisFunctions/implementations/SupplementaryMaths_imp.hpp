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
  * @file SupplementaryMaths.hpp
  *
  * @brief declare additional maths feature not implemented in <cmath> 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef SUPPLEMENTARY_MATHS_IMPLEMENTATIONS 
#define SUPPLEMENTARY_MATHS_IMPLEMENTATIONS 

#ifndef SUPPLEMENTARY_MATHS_DECLARATIONS 
    #include "SupplementaryMaths.hpp" 
#endif 

namespace BasisFunctions {

    template < typename Z >
    Z Factorial ( const Z n ) {

        if ( n < 0 ) {

            throw std::runtime_error (
                "Factorial: does not support negative numbers"
            );
        }

        if ( n == 0 ) {
            return 1;
        }

        return n * Factorial ( n - 1 );

    }

} // BasisFunctions : Factorial 


namespace BasisFunctions {

    template < typename Z, typename C >
    C operator* ( const Z integer, const C complex ) {

        C result;

        result.real ( integer * complex.real() );
        result.imag ( integer * complex.imag() );

        return result;

    }

} // BasisFunctions : integer & complex multiplication 


namespace BasisFunctions {

    template < typename Z >
    Z Binomial ( const Z n, const Z k ) {

        if ( n < 0 || k < 0 ) {

            throw std::runtime_error (
                "Binomial: currently doesn't support negative values"
            );

        }

        if ( n < k ) { return 0; }

        if ( n == k || k == 0 ) {
            return 1;
        }

        return Binomial<Z> ( n - 1, k - 1 ) +
               Binomial<Z> ( n - 1, k     );

    }

} // BasisFunctions : Binomial 


#endif // SUPPLEMENTARY_MATHS_IMPLEMENTATIONS 

