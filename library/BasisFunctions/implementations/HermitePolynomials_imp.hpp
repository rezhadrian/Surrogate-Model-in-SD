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

#ifndef BASIS_FUNCTIONS_DECLARATIONS 
    #include "BasisFunctions.hpp" 
#endif 


namespace BasisFunctions {

    template < typename Z, typename C >
    Vector<C> HermitePolynomials ( 
        const Vector<Z>& indices, const Vector<C>& X, const Z dim 
    ) {

        if ( dim <= 0 ) {

            throw std::runtime_error (
                "HermitePolynomials: dimension must be positive"
            );

        }

        auto nProducts = indices.size() / dim;
        auto nSamples  =       X.size() / dim;

        if ( indices.size() - dim * nProducts != 0 ) {

            throw std::runtime_error (
                "HermitePolynomials: num of indices not multiple of dimension"
            );

        }

        if (       X.size() - dim * nSamples  != 0 ) {

            throw std::runtime_error (
                "HermitePolynomials: num of samples not multiple of dimension"
            );

        }

        Vector<C> result;
        result.reserve ( nProducts * nSamples );

        for ( auto i = 0; i < nSamples ; i++ ) {
        for ( auto j = 0; j < nProducts; j++ ) {

            C unity = 1.0;

            result.push_back ( 
                std::transform_reduce (
                
                    indices.begin() + j * dim,
                    indices.begin() + j * dim + dim,
                          X.begin() + i * dim,

                    unity,

                    // resulting polynomials are multiplied
                    []( const auto N, const auto P ) { return N * P ; },

                    // polynomials with given index and variable
                    []( const auto idx, const auto x ) {

                        return HermitePolynomial<Z,C> ( idx, x ) / 

                               std::sqrt ( Factorial<Z> ( idx ) );

                    }
                ) 
            );

        }
        }

        return result;

    }

} // BasisFunctions : HermitePolynomials 


namespace BasisFunctions {

    template < typename Z, typename C >
    C HermitePolynomial ( const Z index, const C x ) {

        if ( index == 0 ) return 1;
        if ( index == 1 ) return x;

        return (

                        x * HermitePolynomial ( index - 1, x ) -
            ( index - 1 ) * HermitePolynomial ( index - 2, x )

        );     

    }

} // BasisFunctions : HermitePolynomial 


#endif // HERMITE_POLYNOMIALS_IMPLEMENTATIONS 

