/**
  * @file SupplementaryMaths_BF.hpp 
  *
  * @brief declarations of math functions not implemented in cmath 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef SUPPLEMENTARY_MATHS_BF_DECLARATIONS 
#define SUPPLEMENTARY_MATHS_BF_DECLARATIONS 


namespace BasisFunctions {


    template < typename Z > 
    /**
      * Compute n! recursively 
      * @tparam Z a type of non-negative integer e.g. size_t 
      */
    Z Factorial ( const Z n ); 


    template < typename Z, typename C >
    /**
      * Overload multiplication for integer and complex operand 
      *
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam C a class comparable with std::complex 
      */
    C operator* ( const Z integer, const C complex );


    template < typename Z >
    /**
      * Compute binomial coefficient nCk recursively 
      *
      * @tparam Z a type of non-negative integer e.g. size_t 
      */
    Z Binomial ( const Z n, const Z k );


} // BasisFunctions : SupplementaryMaths 


#endif // SUPPLEMENTARY_MATHS_BF_DECLARATIONS 

