/**
  * @file Arithmetic.hpp 
  *
  * @brief declarations of templated arithmetic matrix operations .
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef ARITHMETIC_DECLARATIONS 
#define ARITHMETIC_DECLARATIONS 

#include "LibrariesLoader_LA.hpp" 


template < class Vector, typename T >
/**
  * Multiplication operator between a scalar and a vector.
  *
  * @param c a complex scalar.
  * @param v an object with functions size() and operator[].
  */
Vector operator* ( T c, const Vector& v );


template < class Vector >
/**
  * Element-wise addition between two vectors.
  *
  * @param v1 an object with functions size() and operator[].
  * @param v2 an object with functions size() and operator[].
  *
  * @return the vector (v1 + v2).
  */
Vector operator+ ( const Vector& v1, const Vector& v2 );


template < class Vector >
/**
  * Element-wise subtraction between two vectors.
  *
  * @param v1 an object with functions size() and operator[].
  * @param v2 an object with functions size() and operator[].
  *
  * @return the vector (v1 - v2).
  */
Vector operator- ( const Vector& v1, const Vector& v2 );


#ifndef ARITHMETIC_IMPLEMENTATIONS 
    #include "Arithmetic_imp.hpp" 
#endif 

#endif // ARITHMETIC_DECLARATIONS 

