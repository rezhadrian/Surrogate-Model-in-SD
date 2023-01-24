/**
  * @file SupplementaryMaths_imp.hpp
  *
  * @brief implement additional maths feature not implemented in <cmath> 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef SUPPLEMENTARY_MATHS_IMPLEMENTATIONS 
#define SUPPLEMENTARY_MATHS_IMPLEMENTATIONS 

#ifndef BASIS_FUNCTIONS_DECLARATIONS 
    #include "BasisFunctions.hpp" 
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

