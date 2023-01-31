/**
  * @file Arithmetic_imp.hpp 
  *
  * @brief implementations of templated arithmetic matrix operations .
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef ARITHMETIC_IMPLEMENTATIONS 
#define ARITHMETIC_IMPLEMENTATIONS 

#ifndef ARITHMETIC_DECLARATIONS 
    #include "Arithmetic.hpp" 
#endif 


template < class Vector, typename T>
Vector operator* ( T c, const Vector& v ) {

    Vector result ( v.size() );

    for ( auto i = 0; i < v.size(); i++ ) {

        result[i] = c * v[i];;

    }

    return result;

}

    
template < class Vector > 
Vector operator+ ( const Vector& v1, const Vector& v2 ) {

    if ( v1.size() != v2.size() ) {

        throw std::runtime_error (
            "Vector addition: sizes don't match"
        );

    }

    Vector result ( v1.size() );

    for ( auto i = 0; i < v1.size(); i++ ) {

        result[i] = v1[i] + v2[i];

    }

    return result;

} // operator+ 


template < class Vector > 
Vector operator- ( const Vector& v1, const Vector& v2 ) {

    if ( v1.size() != v2.size() ) {

        throw std::runtime_error (
            "Vector subtraction: sizes don't match"
        );

    }

    Vector result ( v1.size() );

    for ( auto i = 0; i < v1.size(); i++ ) {

        result[i] = v1[i] - v2[i];

    }

    return result;

} // operator- 


#endif // ARITHMETIC_IMPLEMENTATIONS 

