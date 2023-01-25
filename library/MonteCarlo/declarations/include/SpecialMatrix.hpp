/**
  * @file SpecialMatrix.hpp 
  *
  * @brief declare special matrix classes used in Monte Carlo 
  * 
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef SPECIAL_MATRIX_DECLARATIONS 
#define SPECIAL_MATRIX_DECLARATIONS 

#include "LibrariesLoader_MC.hpp" 

template < typename T >
using Vector = std::vector<T>;


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
      * Dense upper triangular matrix 
      * Stores only lower triangular and diagonal elements
      * Stores zeros as well
      * 
      * @tparam Z a type of non-negative integers e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      */
    class DenseUTriangularMatrix {

                Z dim_;
        Vector<R> dat_;

        public: 

            DenseUTriangularMatrix ( const Z dim );

            R  operator() ( const Z i, const Z j ) const;
            R& operator() ( const Z i, const Z j );

            Z dimension () const { return dim_; };

            Vector<R>& data () { return dat_; }

    };


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

            Vector<R>& data () { return dat_; }

    }; // DenseLTriangularMatrix 


} // MonteCarlo : SpecialMatrix 


#endif // SPECIAL_MATRIX_DECLARATIONS 

