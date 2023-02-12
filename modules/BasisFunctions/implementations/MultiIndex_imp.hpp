/**
  * @file MultiIndex_imp.hpp
  *
  * @brief 
  * Implementations of functions to calculate indices of univariate function. 
  * 
  * @anchor _MultiIndex_imp_hpp_ 
  *
  * @author 
  * Rezha Adrian Tanuharja @n 
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
      * @private 
      * 
      * @brief 
      * Compute binomial coefficient nCk recursively. @n 
      * Used instead of boost implementation which return floating number. 
      *
      * @tparam Z a type of non-negative integer e.g. size_t 
      * 
      * @param n size of set to sample from 
      * @param k size of subset 
      * 
      * @return number of distinct k-subset of n-set 
      */
    Z Binomial ( const Z n, const Z k );

} // BasisFunctions : Additional math functions not in cmath 


namespace BasisFunctions {

    template < typename Z >
    /**
      * @private 
      * 
      * @brief 
      * Generate a new set of indices of univariate functions. 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      *
      * @param n      number of positive integers used to generate indices 
      * @param Subset vector to save the new set 
      * @param Active track all previous sets that have been saved.
      */
    void MultiIndexRecursive ( 
        const Z n, Vector<Z>& Subset, Marker& Active 
    );

} // BasisFunctions : MultiIndexRecursive 


namespace BasisFunctions {

    template < typename Z >
    Vector<Z> MultiIndex ( const Z SetSize, const Z iMax ) {

        if ( SetSize <= 0 ) {

            throw std::runtime_error (
                "MultiIndex: dimension must be positive integers"
            );

        }

        if ( iMax < 0 ) {

            throw std::runtime_error (
                "MultiIndex: max order cannot be negative"
            );

        }

        Vector<Z> result ( SetSize, 0 );
        result.reserve ( SetSize * Binomial<Z> (iMax + SetSize, SetSize) );
    
        for ( auto i = 0; i < iMax; i++ ) {

            auto jMax = Binomial<Z> ( SetSize + i, SetSize - 1 );

            Marker Active ( SetSize + i, false );
            std::transform (
                Active.begin(),
                Active.begin() + SetSize - 1,
                Active.begin(),

                []( auto m ) {
                    return m + true;
                }

            );

            for ( auto j = 0; j < jMax; j++ ) {

                MultiIndexRecursive<Z> ( SetSize + i, result, Active );

            }

        }

        return result;

    }

} // BasisFunctions : MultiIndex 


namespace BasisFunctions {

    template < typename Z >
    void MultiIndexRecursive ( 
        const Z n, Vector<Z>& Subset, Marker& Active 
    ) {

        if ( Active.size() != n ) {

            throw std::runtime_error (
                "MultiIndexRecursive: tracker size doesn't match nSet"
            );

        }

        Vector<Z> result;

        Vector<Z> Set ( n );
        std::iota ( Set.begin(), Set.end(), 1 );

        std::next_permutation ( Active.begin(), Active.end() );

        Z m = 0;
        for ( const auto& marker : Active ) {

            if ( marker ) {
                result.push_back ( Set[m] );
            }
            
            m++;

        }

        result.push_back ( n + 1 );

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

