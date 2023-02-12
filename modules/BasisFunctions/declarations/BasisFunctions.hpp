/**
  * @file BasisFunctions.hpp 
  *
  * @brief 
  * Declarations of all functions needed to generate basis functions. 
  * 
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef BASIS_FUNCTIONS_DECLARATIONS 
#define BASIS_FUNCTIONS_DECLARATIONS 

#include "LibrariesLoader_BF.hpp" 

/** 
  * @namespace BasisFunctions 
  * 
  * @brief 
  * Contains all functions needed to generate basis functions. 
  * 
  * @anchor _BasisFunctions_ 
  */
namespace BasisFunctions {

    template < typename T >
    using Vector = std::vector <T>;

    using Marker = std::vector <bool>;


    template < typename Z >
    /**
      * @brief 
      * Generate sets of indices of univariate functions. 
      * Sets with lowest indices sum is the first element in the vector. @n 
      * Implemented in @ref _MultiIndex_imp_hpp_ 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      *
      * @param  SetSize number of indices in a set 
      * @param  iMax    largest allowable index 
      * 
      * @return vector of indices, size of which is a multiple of SetSize 
      */
    Vector<Z> MultiIndex ( const Z SetSize, const Z iMax );


    template < typename Z >
    /**
      * @brief 
      * Remove sets of indices with sum of indices above certain treshold. 
      * Operations are performed in-place, no new vector is generated. @n
      * Implemented in @ref _Truncations_imp_hpp_ 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      *
      * @param Indices vector of indices, size of which is a multiple of SetSize 
      * @param SetSize number of indices in a set 
      * @param MaxSum  largest allowable sum of indices in a set  
      */
    void TotalTruncation ( 
        Vector<Z>& Indices, const Z SetSize, const Z MaxSum 
    );


    template < typename Z, typename R, typename C >
    /**
      * @brief 
      * Evaluate probabilist Hermite polynomials and multiply polynomials
      * from the same set. @n 
      * Implemented in @ref _HermitePolynomials_imp_hpp_ 
      *
      * @tparam Z a type of non-negative integer e.g. size_t
      * @tparam R a type of floating number e.g. double 
      * @tparam C a type of floating complex number e.g. std::complex<float> 
      *
      * @param Indices vector of indices of Hermite polynomial 
      * @param Args    vector of arguments for Hermite polynomial 
      * @param SetSize number of polynomials in a set 
      * 
      * @return vector of products of Hermite polynomials
      */
    Vector<C> HermitePolynomials (
        const Vector<Z>& Indices, const Vector<C>& Args, const Z SetSize  
    );


    template < typename Z, typename R > 
    /**
      * @brief 
      * Compute expected value of products of three Hermite polynomials 
      * for all triplet pairs. @n 
      * Implemented in @ref _TripleHermite_imp_hpp_ 
      * 
      * @tparam Z a type of non-negative integer e.g size_t 
      * @tparam R a type of floating number e.g. double 
      * 
      * @param Indices vector of indices of Hermite polynomial 
      * @param SetSize number of polynomials in a set 
      * @param k       index of the first polynomial in the triplet 
      * 
      * @return vector containing Prod ( E(H_k,H_i,H_j) ) 
      */
    Vector<R> ExpHermiteTriples ( 

        const Vector<Z>& Indices, 
        const Z SetSize, 
        const Z k

    );


} // BasisFunctions 


#ifndef MULTI_INDEX_IMPLEMENTATIONS 
    #include "MultiIndex_imp.hpp" 
#endif 

#ifndef TRUNCATIONS_IMPLEMENTATIONS 
    #include "Truncations_imp.hpp" 
#endif 

#ifndef HERMITE_POLYNOMIALS_IMPLEMENTATIONS 
    #include "HermitePolynomials_imp.hpp" 
#endif 

#ifndef TRIPLE_HERMITE_IMPLEMENTATIONS 
    #include "TripleHermite_imp.hpp" 
#endif 

#endif // BASIS_FUNCTIONS_DECLARATIONS 

