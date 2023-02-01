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


} // BasisFunctions : Additional math functions not in cmath 


namespace BasisFunctions {

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

} // BasisFunctions : HermitePolynomial 


namespace BasisFunctions {

    template < typename Z, typename R, typename C >
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

                               std::sqrt<R> ( Factorial<Z> ( idx ) );

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


#endif // HERMITE_POLYNOMIALS_IMPLEMENTATIONS 

