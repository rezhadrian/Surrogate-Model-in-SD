/**
  * @file MonteCarlo.hpp
  *
  * @brief 
  * Declarations of all  functions required to run computer experiments. 
  *
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef MONTE_CARLO_DECLARATIONS 
#define MONTE_CARLO_DECLARATIONS 

#include "LibrariesLoader_MC.hpp" 

/** 
  * @namespace MonteCarlo 
  * 
  * @brief 
  * Contains all functions needed to run computer experiments. 
  * 
  * @anchor _MonteCarlo_ 
  */
namespace MonteCarlo {

    template < typename T >
    using Vector = std::vector<T>;

    template < typename T >
    using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;


    template < typename Z, typename R >
    /**
      * @brief
      * Perform latin hypercube sampling with given dimension. @n 
      * Implemented in @ref _LatinHypercubeSampling_imp_hpp_ 
      * 
      * @anchor _LHS_ 
      *
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating number e.g. double 
      * 
      * @param nPoints number of sampled points 
      * @param Dim     dimension of each point 
      * 
      * @return vector { Point1, Point2, ... } 
      */
    Vector<R> LHS ( const Z nPoints, const Z Dim );


    template < typename Z, typename R >
    /**
      * @brief 
      * Convert LHS result to sets of standard normal random variables.
      * Conversion is in place, no new vector generated. @n 
      * Implemented in @ref _VariableGeneration_imp_hpp_ 
      *
      * @tparam Z a type of non=negative integer e.g. size_t 
      * @tparam R a type of floating number e.g. double 
      * 
      * @param LHSResult vector output of @ref _LHS_ 
      */
    void ConvertLHStoStdNorm ( Vector<R>& LHSResult );


    template < typename R, typename C >
    /**
      * @brief 
      * Combine uncorrelated RVs to generate correlated RVs. @n 
      * Implemented in @ref _VariableGeneration_imp_hpp_ 
      *
      * @tparam R a type of floating number e.g. double 
      * @tparam C a type of floating complex number e.g. std::complex<float> 
      * 
      * @param Correl correlation matrix 
      * @param RVs    vector of uncorrelated random variables 
      * @param Dim    number of random variables that will be correlated 
      * 
      * @return vector of correlated random variables casted as complex 
      */
    Vector<C> CombineRVs ( 

        const MatrixXT<R>& Correl, 
        Vector<R>& RVs,
        const size_t Dim

    );


    template < typename Z, typename R >
    /**
      * @brief 
      * Generate RVs with a given distribution from standard normal RVs. @n 
      * Implemented in @ref _VariableGeneration_imp_hpp_ 
      * 
      * @tparam Z a type of non-negative integer e.g. size_t 
      * @tparam R a type of floating number e.g. double 
      * 
      * @param StdNormRVs vector of standard normal random variables 
      * @param Correl     correlation matrix 
      * @param ICDFs      vector of inverse cummulative distribution functions 
      * @param Dim        number of distributions 
      * 
      * @return vector of new randow variables 
      */
    Vector<R> GenerateRVs (

        Vector<R>& StdNormRVs, 
        const MatrixXT<R>& Correl, 
        const Vector< std::function<R(R)> >& ICDFs,
        const Z Dim

    );


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

