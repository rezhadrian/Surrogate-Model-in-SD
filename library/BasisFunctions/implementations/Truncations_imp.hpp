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
  * @file Truncations_imp.hpp
  *
  * @brief implement functions to truncate indices of orthonormal polynomials.
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef TRUNCATIONS_IMPLEMENTATIONS 
#define TRUNCATIONS_IMPLEMENTATIONS 

#ifndef TRUNCATIONS_DECLARATIONS 
    #include "Truncations.hpp" 
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

