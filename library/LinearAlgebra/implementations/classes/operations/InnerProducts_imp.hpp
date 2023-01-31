/**
  * @file InnerProducts_imp.hpp 
  *
  * @brief implementations of templated vector inner products 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezhadr@outlook.com 
  */

#ifndef INNER_PRODUCTS_IMPLEMENTATIONS 
#define INNER_PRODUCTS_IMPLEMENTATIONS 

#ifndef INNER_PRODUCTS_DECLARATIONS 
    #include "InnerProducts.hpp" 
#endif 


namespace linalg {


    template < class VectorType1, class VectorType2, typename T >
    T Euclidean ( const VectorType1& v1, const VectorType2& v2 ) {

        if ( v1.size() != v2.size() ) {

            throw std::runtime_error (
                "Euclidean InnerProduct: sizes don't match"
            );

        }

        T result = 0.0;

        for ( auto i = 0; i < v1.size(); i++ ) {

            result += std::conj(v1[i]) * v2[i];

        }

        return result;

    }


} // linalg 


#endif // INNER_PRODUCTS_IMPLEMENTATIONS 

