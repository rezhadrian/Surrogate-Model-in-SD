/**
  * @file TripleHermite_imp.hpp
  *
  * @brief 
  * Implementations of functions to calculate exp value of Hermite triples 
  * 
  * @anchor _TripleHermite_imp_hpp_ 
  *
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef TRIPLE_HERMITE_IMPLEMENTATIONS 
#define TRIPLE_HERMITE_IMPLEMENTATIONS 

#ifndef BASIS_FUNCTIONS_DECLARATIONS 
    #include "BasisFunctions.hpp" 
#endif 

#include <boost/math/special_functions/factorials.hpp> 


namespace BasisFunctions {

    template < typename Z, typename R > 
    /**
      * @private 
      * 
      * @brief 
      * Compute expected value of products of three Hermite polynomials. 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating number e.g. double 
      * 
      * @param i index of the first polynomial 
      * @param j index of the second polynomial 
      * @param k index of the third polynomial 
      * 
      * @return E( H_i, H_j, H_k ) 
      */
    R EHermiteTriple ( const Z i, const Z j, const Z k );

} // BasisFunctions : EHermiteTriple 


namespace BasisFunctions {

    template < typename Z, typename R >
    Vector<R> ExpHermiteTriples ( 

        const Vector<Z>& indices, 
        const Z dim, 
        const Z k

    ) {

        auto nPoints = indices.size() / dim;

        Vector<R> result;
        result.reserve ( nPoints * nPoints );


        for ( auto i = 0; i < nPoints; i++ ) {
        for ( auto j = 0; j < nPoints; j++ ) {

            result.push_back (

                std::transform_reduce (

                    indices.begin() + i * dim, 
                    indices.begin() + i * dim + dim, 
                    indices.begin() + j * dim, 

                    1.0, 

                    // Expected values are multiplied point-wise 
                    []( const auto E1, const auto E2 ) { return E1 * E2; }, 

                    // Expected value of triple Hermite product 
                    [k]( const auto m, const auto n ) {

                        return EHermiteTriple <Z,R> ( m, n, k );

                    }

                )

            );

        }
        }

        return result; 

    }

} // BasisFunctions : EHermiteTriples 


namespace BasisFunctions {

    template < typename Z, typename R > 
    R EHermiteTriple ( const Z i, const Z j, const Z k ) {

        auto s = i + j + k;

        if ( s % 2 != 0 )                  return 0.0;
        if ( s < 2 * std::max( {i,j,k} ) ) return 0.0;

        return 

            std::sqrt ( 

                boost::math::factorial<R> (i) * 
                boost::math::factorial<R> (j) * 
                boost::math::factorial<R> (k) 

            ) / (

                boost::math::factorial<R> ( s/2 - i ) * 
                boost::math::factorial<R> ( s/2 - j ) * 
                boost::math::factorial<R> ( s/2 - k ) 

            );

    }

} // BasisFunctions : EHermiteTriple  


#endif // TRIPLE_HERMITE_IMPLEMENTATIONS 

