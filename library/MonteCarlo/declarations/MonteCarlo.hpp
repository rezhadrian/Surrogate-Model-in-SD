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
  * @file MonteCarlo.hpp
  *
  * @brief declare functions required to run computer experiments 
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef MONTE_CARLO_DECLARATIONS 
#define MONTE_CARLO_DECLARATIONS 

#include "LibrariesLoader_MC.hpp" 


template < typename T >
using Vector = std::vector <T>;


// Implemented in SupplementaryMaths_imp_MC.hpp 

namespace MonteCarlo {


    template < typename Z, typename R >
    /**
      * Raise number to given power recursively 
      *
      * @tparam Z a type of non-negative integers e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      */
    R Power ( const R number, const Z power );


} // MonteCarlo : SupplementaryMaths 


// Implemented in LatinHypercubeSampling_imp.hpp 

namespace MonteCarlo {


    template < typename Z, typename R,
               
               class RandomDevice,
               class RandomEngine,
               class Distribution,
               class Shuffler 

    >
    /**
      * Perform latin hypercube sampling with given sizes 
      *
      * @tparam Z a type of non=negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      * 
      * @tparam RandomDevice non-deterministic random number generator
      *         e.g. std::random_device 
      * @tparam RandomEngine pseudo-random number generator between 0 and 1
      *         e.g. std::default_random_engine 
      * @tparam Distribution distribution within each interval 
      *         e.g. std::uniform_real_distribution<double> 
      * @tparam Shuffler pseudo random shuffler 
      *         e.g. std::mt19937
      * 
      * @param nInterval number of intervals in [0,1] 
      * @param nSample   number of sampled [0,1] 
      * @return vector { Sample1, Sample2, ... } 
      */
    Vector<R> LHS ( const Z nInterval, const Z nSample );


} // MonteCarlo : LatinHypercubeSampling 


// Implemented in InverseCDF_imp.hpp 

namespace MonteCarlo {


    template < typename Z, typename R >
    /**
      * Apply Inverse CDF function to a vector of uniformly distributed RVs
      * Results overwrite inputs 
      *
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      *
      * @param UniformSamples vector with elements in (0,1) 
      */
    void InverseCDFs ( Vector<R>& UniformSamples );


    template < typename Z, typename R >
    /**
      * Compute inverse of standard normal CDF for a given quantile 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      * 
      * @return z such that CDF(z) = quantile 
      */
    R InvStdNormCDF ( const R quantile );


    template < typename R >
    /**
      * Compute inverse of complement error function 
      *
      * @tparam R a type of floating point e.g. double 
      * @return x such that erfc(z) = x  
      */
    R InvErfc ( const R p );


    template < typename R >
    /**
      * Compute inverse of log normal CDF for a given quantile 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      *
      * @param mu mean of distribution 
      * @param sig standard deviation of distribution 
      * 
      * @return z such that CDF(z) = quantile 
      */
    R InvLogNormCDF ( const R quantile, const R mu, const R sig );


} // MonteCarlo : InverseCDF 


// Implemented in SpecialMatrix_imp.hpp 

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


// Implemented in CholeskyDecomposition_imp.hpp 

namespace MonteCarlo {


    template < class LTriangularMatrix, class SymMatrix >
    /**
      * Cholesky decomposition of SPD matrix e.d. covariance matrix 
      *
      * @tparam LTriangularMatrix lower triangular matrix class 
      * @tparam SymMatrix symmetric matrix class 
      * 
      * @return lower triangular part of the decomposition 
      */
    LTriangularMatrix CholeskyDecompose ( const SymMatrix& A );


} // MonteCarlo : CholeskyDecomposition  


#ifndef SUPPLEMENTARY_MATHS_MC_IMPLEMENTATIONS 
    #include "SupplementaryMaths_imp_MC.hpp" 
#endif 

#ifndef LATIN_HYPERCUBE_SAMPLING_IMPLEMENTATIONS 
    #include "LatinHypercubeSampling_imp.hpp" 
#endif 

#ifndef INVERSE_CDF_IMPLEMENTATIONS 
    #include "InverseCDF_imp.hpp" 
#endif 

#ifndef SPECIAL_MATRIX_IMPLEMENTATIONS 
    #include "SpecialMatrix_imp.hpp" 
#endif 

#ifndef CHOLESKY_DECOMPOSITION_IMPLEMENTATIONS 
    #include "CholeskyDecomposition_imp.hpp" 
#endif 

#endif // MONTE_CARLO_DECLARATIONS 

