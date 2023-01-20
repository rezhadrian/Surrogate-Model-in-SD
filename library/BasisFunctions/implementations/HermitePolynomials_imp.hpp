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
  * @file HermitePolynomials_imp.hpp
  *
  * @brief implement functions to calculate probabilist Hermite polynomials.
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef HERMITE_POLYNOMIALS_IMPLEMENTATIONS 
#define HERMITE_POLYNOMIALS_IMPLEMENTATIONS 

#ifndef HERMITE_POLYNOMIALS_DECLARATIONS 
    #include "HermitePolynomials.hpp" 
#endif 


namespace BasisFunctions {

    template < typename Z, typename C >
    Vector<C> HermitePolynomials ( Vector<Z>& indices, Vector<C>& X, Z dim ) {

        auto nProducts = indices.size() / dim;

        Vector<C> result;
        result.reserve ( nProducts );

        Vector<C> tupple ( dim );

        for ( auto i = 0; i < nProducts; i++ ) {

            std::transform (

                indices.begin() + i * dim,
                indices.begin() + i * dim + dim,
                      X.begin() + i * dim,

                 tupple.begin(),

                 []( auto n, auto x ) {
                    return HermitePolynomial<Z,C> ( n, x ) /
                           std::sqrt ( Factorial<Z> ( n ) );
                 }

            );

            C unity = 1.0;

            result.push_back (
                std::accumulate (
                    tupple.begin(), tupple.end(),
                    unity,
                    std::multiplies<C>()
                )
            );

        }

        return result;

    }

} // BasisFunctions : HermitePolynomials 


namespace BasisFunctions {

    template < typename Z, typename C >
    C HermitePolynomial ( Z index, C x ) {

        if ( index == 0 ) return 1;
        if ( index == 1 ) return x;

        return (

                                    x * HermitePolynomial ( index - 1, x ) -
            Multiply<Z,C> ( index - 1,  HermitePolynomial ( index - 2, x ) )

        );     

    }

} // BasisFunctions : HermitePolynomial 


namespace BasisFunctions {

    template < typename Z, typename C > 
    C Pow ( Z power, C x ) {

        C result = 1.0;

        for ( auto i = 0; i < power; i++ ) {
            result *= x;
        }

        return result;

    }

} // BasisFunctions : Pow 


namespace BasisFunctions {

    template < typename Z >
    Z Factorial ( Z n ) {

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
    C Multiply ( Z index, C x ) {

        C result;

        result.real ( index * x.real() );
        result.imag ( index * x.imag() );

        return result;

    }

} // BasisFunctions : Multiply 


#endif // HERMITE_POLYNOMIALS_IMPLEMENTATIONS 

