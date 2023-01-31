/**
  * @file InnerProducts.hpp 
  *
  * @brief declarations of templated vector inner products 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef INNER_PRODUCTS_DECLARATIONS 
#define INNER_PRODUCTS_DECLARATIONS 

#include "LibrariesLoader_LA.hpp" 

namespace linalg {


    template < class VectorType1, class VectorType2, typename T >
    /**
      * Compute Euclidean inner product between two vectors 
      * 
      * @param v1 an object with functions size() and operator[].
      * @param v2 an object with functions size() and operator[].
      * 
      * @return <v1,v2> 
      */
    T Euclidean ( const VectorType1& v1, const VectorType2& v2 );


} // linalg 

#ifndef INNER_PRODUCTS_IMPLEMENTATIONS 
    #include "InnerProducts_imp.hpp" 
#endif 

#endif // INNER_PRODUCTS_DECLARATIONS 

