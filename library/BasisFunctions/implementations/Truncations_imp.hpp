/**
  * @file Truncations_imp.hpp
  *
  * @brief 
  * Implementations of functions to truncate indices of univariate function.
  *
  * @anchor _Truncations_imp_hpp_ 
  * 
  * @author 
  * Rezha Adrian Tanuharja @n 
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
        Vector<Z>& Indices, const Z SetSize, const Z MaxSum 
    ) {

        auto nSet = Indices.size();

        if ( SetSize * (nSet / SetSize) - nSet != 0 ) {

            throw std::runtime_error (
                "Total truncation: Indices size not a multiple of dimension"
            );

        }

        auto i = 0;
        while ( true ) {

            if ( i >= nSet ) break;

            if ( std::accumulate (

                    Indices.begin() + i,
                    Indices.begin() + i + SetSize,
                    0

                 ) <= MaxSum 

            ) { i += SetSize; continue; }

            Indices.erase (
                Indices.begin() + i,
                Indices.begin() + i + SetSize
            );

            nSet -= SetSize;

        }

    } 

} // BasisFunctions : TotalTruncation 


#endif // TRUNCATIONS_IMPLEMENTATIONS 

