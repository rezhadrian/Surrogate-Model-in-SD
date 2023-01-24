/**
  * @file Truncations_imp.hpp
  *
  * @brief implement functions to truncate indices of orthonormal polynomials.
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef TRUNCATIONS_IMPLEMENTATIONS 
#define TRUNCATIONS_IMPLEMENTATIONS 

#ifndef BASIS_FUNCTIONS_DECLARATIONS 
    #include "BasisFunctions.hpp" 
#endif 


namespace BasisFunctions {

    template < typename Z >
    void TotalTruncation ( 
        Vector<Z>& MultiIndex, const Z dim, const Z SumPMax 
    ) {

        auto nTupple = MultiIndex.size();

        if ( dim * (nTupple / dim) - nTupple != 0 ) {

            throw std::runtime_error (
                "Total truncation: Indices size not a multiple of dimension"
            );

        }

        auto i = 0;
        while ( true ) {

            if ( i >= nTupple ) break;

            if ( std::accumulate (

                    MultiIndex.begin() + i,
                    MultiIndex.begin() + i + dim,
                    0

                 ) <= SumPMax 

            ) { i += dim; continue; }

            MultiIndex.erase (
                MultiIndex.begin() + i,
                MultiIndex.begin() + i + dim
            );

            nTupple -= dim;

        }

    } 

} // BasisFunctions : TotalTruncation 


#endif // TRUNCATIONS_IMPLEMENTATIONS 

