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
  * @file MultiIndex_imp.hpp
  *
  * @brief implement functions to calculate indices of orthonormal polynomials.
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef MULTI_INDEX_IMPLEMENTATIONS 
#define MULTI_INDEX_IMPLEMENTATIONS 

#ifndef MULTI_INDEX_DECLARATIONS 
    #include "MultiIndex.hpp" 
#endif 


namespace BasisFunctions {

    template < typename Z >
    Vector<Z> MultiIndex ( Z dim, Z pMax ) {

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
    void MultiIndexRecursive ( Z nSet, Vector<Z>& Subset, Marker& Active ) {

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
    Z Binomial ( Z DimSet, Z DimSubset ) {

        if ( DimSet < 0 || DimSubset < 0 ) {

            throw std::runtime_error (
                "Binomial: currently doesn't support negative values"
            );

        }

        if ( DimSet < DimSubset ) { return 0; }

        if ( DimSet == DimSubset || DimSubset == 0 ) {
            return 1;
        }

        return Binomial<Z> ( DimSet - 1, DimSubset - 1 ) +
               Binomial<Z> ( DimSet - 1, DimSubset     );

    }

} // BasisFunctions : Binomial 


#endif // MULTI_INDEX_IMPLEMENTATIONS 

