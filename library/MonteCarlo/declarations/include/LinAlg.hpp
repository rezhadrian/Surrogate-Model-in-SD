/**
  * @file LinAlg.hpp 
  *
  * @brief declare linear algebra functionalities used in MonteCarlo 
  * 
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef LINALG_DECLARATIONS 
#define LINALG_DECLARATIONS 

#include <vector> 
#include <cmath> 

template < typename T >
using Vector = std::vector<T>;


// Special matrix classes 

namespace MonteCarlo {


    template < typename Z, typename R >
    /**
      * Dense symmetric matrix 
      * Stores only lower triangular and diagonal elements
      * Stores zeros as well
      * 
      * @tparam Z a type of non-negative integers e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      */
    class DenseSymMatrix {

                Z dim_;
        Vector<R> dat_;

        public: 

            DenseSymMatrix ( const Z dimension );

            R  operator() ( const Z i, const Z j ) const;
            R& operator() ( const Z i, const Z j );

            Z dimension () const { return dim_; };

            Vector<R>& data () { return dat_; }


    }; // DenseSymMatrix 


    template < typename Z, typename R >
    /**
      * Dense lower triangular matrix 
      * Stores only lower triangular and diagonal elements
      * Stores zeros as well
      * 
      * @tparam Z a type of non-negative integers e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      */
    class DenseLTriangularMatrix {

                Z dim_;
        Vector<R> dat_;

        public: 

            DenseLTriangularMatrix ( const Z dimension );

            R  operator() ( const Z i, const Z j ) const;
            R& operator() ( const Z i, const Z j );

            Z dimension () const { return dim_; };

            const Vector<R>& data () { return dat_; }

    }; // DenseLTriangularMatrix 


} // MonteCarlo : Special Matrix Classes 


// Cholesky Decomposition 

namespace MonteCarlo {


    template < class LTriangularMatrix, class SymMatrix >
    /**
      * Cholesky decomposition of SPD matrix e.g. correlation matrix 
      *
      * @tparam LTriangularMatrix lower triangular matrix class 
      * @tparam SymMatrix symmetric matrix class 
      * 
      * @return lower triangular part of the decomposition 
      */
    LTriangularMatrix CholeskyDecompose ( const SymMatrix& A );


} // MonteCarlo : CholeskyDecomposition  


#ifndef LINALG_IMPLEMENTATIONS 
    #include "LinAlg_imp.hpp" 
#endif 

#endif // LINALG_DECLARATIONS 

