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

namespace MonteCarlo {

    template < typename T >
    using Vector = std::vector<T>;

}


// Implemented in LatinHypercubeSampling_imp.hpp 

namespace MonteCarlo {


    template < typename Z, typename R >
    /**
      * Perform latin hypercube sampling with given sizes 
      *
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      * 
      * @param nInterval number of intervals in [0,1] 
      * @param nSample   number of sampled [0,1] 
      * @return vector { Sample1, Sample2, ... } 
      */
    Vector<R> LHS ( const Z nInterval, const Z nSample );


} // MonteCarlo : LatinHypercubeSampling 


// Implemented in VariableGeneration_imp.hpp 

namespace MonteCarlo {

    template < typename T >
    using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;


    template < typename Z, typename R >
    /**
      * Convert LHS result to sets of standard normal random variables 
      *
      * @tparam Z a type of non=negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      */
    void ConvertLHStoStdNorm ( Vector<R>& LHSResult );


    template < typename R, typename C >
    /**
      * Combine uncorrelated RVs to generate correlated RVs 
      *
      * @tparam R a type of floating point e.g. double 
      * @tparam LTMatrix lower triangular matrix class 
      * 
      * @param Correl correlation matrix 
      */
    Vector<C> CombineRVs ( 

        const MatrixXT<R>& Correl, 
        Vector<R>& RVs,
        const size_t dim

    );


    template < typename Z, typename R >
    /**
      * Convert standard normal RVs to correlated RVs of other distribution 
      * 
      * @tparam R a type of floating point e.g. double 
      * @tparam LTMatrix a lower triangular matrix 
      * 
      * @param Correl correlation matrix 
      */
    Vector<R> GenerateRVs (

        Vector<R>& StdNormRVs, 
        const MatrixXT<R>& Correl, 
        const Vector< std::function<R(R)> >& ICDFs,
        const Z dim

    );


} // MonteCarlo : ConvertLHS 


namespace MonteCarlo {


    template < typename R, typename C >
    /**
      * Example model to evaluate FRF for a given E and Mu 
      * 
      * @tparam R a type of floating point e.g. double 
      * @tparam C a type compatible with std::complex<R> 
      * 
      * @param omega excitation angular velocity 
      * @param params a vector { E, Mu } 
      */
    C EvaluateFRF ( const R omega, const Vector<R>& params );


    template < typename Z, typename R, typename C, class Function >
    /**
      * Evaluate model result for a set of parameters set 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating point e.g. double 
      * @tparam C a type compatible with std::complex<R> 
      * @tparam Function compatible with C(const R, const Vector<R>&)
      *
      * @param dim number of variable argument for the model
      */
    Vector<C> EvaluateModel ( 

        const R omega, 
        const Vector<R>& RVs, 
        const Function& Model, 
        const Z dim 

    );


} // MonteCarlo : EvaluateFRF 


#ifndef LATIN_HYPERCUBE_SAMPLING_IMPLEMENTATIONS 
    #include "LatinHypercubeSampling_imp.hpp" 
#endif 

#ifndef VARIABLE_GENERATION_IMPLEMENTATIONS 
    #include "VariableGeneration_imp.hpp" 
#endif 

#ifndef ANALYTICAL_MODEL_IMPLEMENTATIONS 
    #include "AnalyticalModel_imp.hpp" 
#endif 

#endif // MONTE_CARLO_DECLARATIONS 

