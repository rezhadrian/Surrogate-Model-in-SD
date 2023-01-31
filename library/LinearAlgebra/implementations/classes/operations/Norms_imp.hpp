/**
  * @file Norms_imp.hpp 
  *
  * @brief implementations of templated vector norms 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef NORMS_IMPLEMENTATIONS 
#define NORMS_IMPLEMENTATIONS 

#ifndef NORMS_DECLARATIONS 
    #include "Norms.hpp" 
#endif 

namespace linalg {


    template < class Vector, typename T >
    T Euclidean ( const Vector& v ) {

        T result = 0.0;

        for ( auto i = 0; i < v.size(); i++ ) {

            result += std::abs(v[i]) * std::abs(v[i]);

        }

        return std::sqrt ( result );

    } // Euclidean 


    template < class Vector, typename T >
    T Max ( const Vector& v ) {

        T result = 0.0;

        for ( auto i = 0; i < v.size(); i++ ) {

            if ( std::abs(v[i]) > result ) {
                result = std::abs(v[i]);
            }

        }

        return result;

    } // Max 


} // linalg 

#endif // NORMS_IMPLEMENTATIONS 

