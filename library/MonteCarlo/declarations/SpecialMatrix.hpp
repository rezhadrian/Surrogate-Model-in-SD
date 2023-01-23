/**
  * This file is part of a Surrogate Model library for structural dynamics. 
  *
  * The library is free: you can distribute it and/or modify it.
  * The library is distributed in the hope that it can be useful and helpful,
  * particularly for students and self - learning individuals.
  * 
  * The library comes WITHOUT ANY WARRANTY: without even the implied warranty 
  * of MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE.
  *
  * The library is licensed under the GNU General Public License v3.0
  *
  *
  * The library is based on works by:
  * - Felix Schneider
  * - Iason Papaioannou
  * 
  * Chair of Structural Mechanics 
  * Engineering Risk Analysis Group 
  *
  * Technische Universitaet Muenchen
  *
  * www.cee.ed.tum.de/bm
  * www.cee.ed.tum.de/era 
  */ 

/**
  * @file SpecialMatrix.hpp
  *
  * @brief declare special matrix class for Monte Carlo  
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef SPECIAL_MATRIX_DECLARATIONS 
#define SPECIAL_MATRIX_DECLARATIONS

#include "LibrariesLoader_MC.hpp" 

namespace MonteCarlo {

    template < typename T >
    using Vector = std::vector <T>;
    

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
      * Dense triangular matrix 
      * Stores only lower triangular and diagonal elements
      * Stores zeros as well
      * 
      * @tparam Z a type of non-negative integers e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      */
    class DenseTriangularMatrix {

                Z dim_;
        Vector<R> dat_;

        public: 

            DenseTriangularMatrix ( const Z dimension );

            R  operator() ( const Z i, const Z j ) const;
            R& operator() ( const Z i, const Z j );

            Z dimension () const { return dim_; };

            Vector<R>& data () { return dat_; }


    }; // DenseTriangularMatrix 


} // MonteCarlo 

#ifndef SPECIAL_MATRIX_IMPLEMENTATIONS 
    #include "SpecialMatrix_imp.hpp" 
#endif 

#endif // SPECIAL_MATRIX_DECLARATIONS 

