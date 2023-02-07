/**
  * @file TripleHermite_imp.hpp
  *
  * @brief implement expectations of products of three Hermite polynomials 
  *
  * @author Rezha Adrian Tanuharja
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
    R ETripleHermite ( const Z i, const Z j, const Z k ) {

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


}



#endif // TRIPLE_HERMITE_IMPLEMENTATIONS 

