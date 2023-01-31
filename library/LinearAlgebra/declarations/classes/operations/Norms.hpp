/**
  * @file Norms.hpp 
  *
  * @brief declarations of templated vector norms 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezhadr@outlook.com 
  */

#ifndef NORMS_DECLARATIONS 
#define NORMS_DECLARATIONS 

#include "LibrariesLoader_LA.hpp" 

namespace linalg {


    template < class Vector, typename T >
    /**
      * Norm induced by dot product in Euclidean space.
      *
      * @param v an object with functions size() and operator[].
      * @return sqrt of dot product between v and itself.
      */
    T Euclidean ( const Vector& v );

    template < class Vector, typename T >
    /**
      * Also known as infinite norm.
      *
      * @param v an object with functions size() and operator[].
      * @return largest absolute value of v's elements.
      */
    T Max ( const Vector& v );

    
} // linalg 


#ifndef NORMS_IMPLEMENTATIONS 
    #include "Norms_imp.hpp" 
#endif 

#endif // NORMS_DECLARATIONS 

