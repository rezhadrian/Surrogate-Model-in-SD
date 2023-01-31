/**
  * @file InnerProducts.hpp 
  *
  * @brief declarations of templated vector inner products 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezhadr@outlook.com 
  */

#ifndef INNER_PRODUCTS_DECLARATIONS 
#define INNER_PRODUCTS_DECLARATIONS 

#include "LibrariesLoader_LA.hpp" 

namespace linalg {


    template < class VectorType1, class VectorType2, typename T >
    T Euclidean ( const VectorType1& v1, const VectorType2& v2 );


} // linalg 

#ifndef INNER_PRODUCTS_IMPLEMENTATIONS 
    #include "InnerProducts_imp.hpp" 
#endif 

#endif // INNER_PRODUCTS_DECLARATIONS 

