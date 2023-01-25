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

// Additional math function not implemented in cmath 
#ifndef SUPPLEMENTARY_MATHS_MC_DECLARATIONS 
    #include "SupplementaryMaths_MC.hpp" 
#endif 

// Special matrix to apply inverse theorem to sets of RVs
#ifndef SPECIAL_MATRIX_DECLARATIONS 
    #include "SpecialMatrix.hpp" 
#endif 


template < typename T >
using Vector = std::vector <T>;


// Implemented in LatinHypercubeSampling_imp.hpp 

namespace MonteCarlo {


    template < typename Z, typename R,
               
        class RandomDevice,   // e.g. std::random_device 
        class RandomEngine,   // e.g. std::default_random_engine 
        class Distribution,   // e.g. std::uniform_real_distribution 
        class Shuffler        // e.g. std::mt19937 

    >
    /**
      * Perform latin hypercube sampling with given sizes 
      *
      * @tparam Z a type of non=negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      * 
      * @param nInterval number of intervals in [0,1] 
      * @param nSample   number of sampled [0,1] 
      * @return vector { Sample1, Sample2, ... } 
      */
    Vector<R> LHS ( const Z nInterval, const Z nSample );


} // MonteCarlo : LatinHypercubeSampling 


// Implemented in InverseCDF_imp.hpp 

namespace MonteCarlo {


    template < typename Z, typename R, typename Function >
    /**
      * Apply given Inverse CDF function to a vector of RVs
      * Results overwrites inputs 
      *
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      * @tparam Function takes one floating point and return another 
      */
    void InverseCDFs ( Vector<R>& RandomVariables, const Function ICDF );


} // MonteCarlo : InverseCDF 


// Implemented in CholeskyDecomposition_imp.hpp 

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

