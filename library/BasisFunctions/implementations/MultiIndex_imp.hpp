/**
  * @file MultiIndex_imp.hpp
  *
  * @brief implement functions to calculate indices of orthonormal polynomials.
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef MULTI_INDEX_IMPLEMENTATIONS 
#define MULTI_INDEX_IMPLEMENTATIONS 

#ifndef BASIS_FUNCTIONS_DECLARATIONS 
    #include "BasisFunctions.hpp" 
#endif 


namespace BasisFunctions {

    template < typename Z >
    /**
      * Compute binomial coefficient nCk recursively 
      *
      * @tparam Z a type of non-negative integer e.g. size_t 
      */
    Z Binomial ( const Z n, const Z k );

} // BasisFunctions : Additional math functions not in cmath 


namespace BasisFunctions {

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

} // BasisFunctions : MultiIndexRecursive 


namespace BasisFunctions {

    template < typename Z >
    Vector<Z> MultiIndex ( const Z dim, const Z pMax ) {

        if ( dim <= 0 ) {

            throw std::runtime_error (
                "MultiIndex: dimension must be positive integers"
            );

        }

        if ( pMax < 0 ) {

            throw std::runtime_error (
                "MultiIndex: max order cannot be negative"
            );

        }

        Vector<Z> result ( dim, 0 );
        result.reserve ( dim * Binomial<Z> (pMax + dim, dim) );
    
        for ( auto i = 0; i < pMax; i++ ) {

            auto jMax = Binomial<Z> ( dim + i, dim - 1 );

            Marker Active ( dim + i, false );
            std::transform (
                Active.begin(),
                Active.begin() + dim - 1,
                Active.begin(),

                []( auto m ) {
                    return m + true;
                }

            );

            for ( auto j = 0; j < jMax; j++ ) {

                MultiIndexRecursive<Z> ( dim + i, result, Active );

            }

        }

        return result;

    }

} // BasisFunctions : MultiIndex 


namespace BasisFunctions {

    template < typename Z >
    void MultiIndexRecursive ( 
        const Z nSet, Vector<Z>& Subset, Marker& Active 
    ) {

        if ( Active.size() != nSet ) {

            throw std::runtime_error (
                "MultiIndexRecursive: tracker size doesn't match nSet"
            );

        }

        Vector<Z> result;

        Vector<Z> Set ( nSet );
        std::iota ( Set.begin(), Set.end(), 1 );

        std::next_permutation ( Active.begin(), Active.end() );

        Z m = 0;
        for ( const auto& marker : Active ) {

            if ( marker ) {
                result.push_back ( Set[m] );
            }
            
            m++;

        }

        result.push_back ( nSet + 1 );

        std::adjacent_difference ( 
            result.begin(), result.end(),
            result.begin()
        );

        std::transform (
            result.begin(), result.end(),
            result.begin(),
            []( auto m ) { return m -1; }
        );

        Subset.insert ( Subset.end(), result.begin(), result.end() );

    }

} // BasisFunctions : MultiIndexRecursive 


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


#endif // MULTI_INDEX_IMPLEMENTATIONS 

