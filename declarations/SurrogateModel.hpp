/**
  * @file SurrogateModel.hpp 
  *
  * @brief 
  * Declarations of collection of surrogate models 
  * 
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef SURROGATE_MODEL_DECLARATIONS 
#define SURROGATE_MODEL_DECLARATIONS 

#include "BasisFunctions.hpp" 
#include "MonteCarlo.hpp" 
#include "AnalyticalModel.hpp"
#include "LibrariesLoader_SM.hpp" 

/**
  * @namespace Surrogate 
  * 
  * @brief 
  * Contains a collection of surrogate structural dynamic model 
  * 
  * @anchor _Surrogate_ 
  */
namespace Surrogate {

    typedef size_t          Z;
    typedef double          R;
    typedef std::complex<R> C; 

    typedef std::vector<R> VectorR; 
    typedef std::vector<C> VectorC; 
    typedef std::vector<Z> VectorZ; 

    typedef Analytical::Model<Z,R,C> AnalyticalModel;

    typedef Eigen::Matrix<C, Eigen::Dynamic, Eigen::Dynamic> MatrixXC;
    typedef Eigen::Vector<C, Eigen::Dynamic> VectorXC;


    /**
      * @class Model 
      *
      * @brief 
      * Abstract class of surrogate structural dynamic model. 
      * 
      * @anchor _Surrogate_Model_ 
      */
    class Model {

        public: 

        virtual Z Dim () const = 0;

        /**
          * @brief 
          * Evaluate surrogate model output for a given parameter. @n 
          * Pure virtual function. 
          * 
          * @param X parameters to evaluate surrogate model 
          * @return approximate output of analytical model 
          */
        virtual VectorC ComputeResponse ( const VectorC& X ) const = 0;

        /**
          * @brief 
          * Train surrogate model with given training data set. @n 
          * Pure virtual function. 
          * 
          * @param TrainSet training data set, usually from LHS 
          */
        virtual void Train ( const VectorC& Load, const VectorR& TrainSet ) = 0;

    }; // Model 


    /**
      * @class NonIntrusivePCE  
      *
      * @brief 
      * Non-intrusive surrogate model based on collocation method. @n 
      * Implemented in @ref _NonIntrusivePCE_cpp_ 
      */
    class NonIntrusivePCE final : public Model {

        const AnalyticalModel* SDModel_;
        VectorZ Indices_;
        VectorC Coeffs_;

        R omega_;
        Z Dim_;

        public: 

        Z Dim () const override { return Dim_; }

        /**
          * @brief 
          * Create PCE model for given SD model and angular velocity 
          * 
          * @param SDModel address of analytical structural dynamic model 
          * @param omega   angular velocity of harmonic load 
          */
        NonIntrusivePCE ( const AnalyticalModel* SDModel, const R omega );

        /**
          * @brief 
          * Set indices of Hermite polynomials to create basis functions 
          * 
          * @param iMax   largest allowable individual index 
          * @param MaxSum largest allowable sum of indices in a set 
          */
        void SetIndices ( Z iMax, Z MaxSum );

        /**
          * @brief 
          * Approximate response of analytical model to given RVs 
          * 
          * @param X random variables added to model's springs 
          * 
          * @return harmonic displacement vector 
          */
        VectorC ComputeResponse ( const VectorC& X ) const override;

        /**
          * @brief 
          * Compute coefficients of basis functions using given training set 
          * 
          * @param Load     harmonic load vector 
          * @param TrainSet vector of random variables 
          */
        void Train ( const VectorC& Load, const VectorR& TrainSet ) override; 

        /**
          * @private 
          */
        void PrintCoeffs () const;
        
    };

} // Surrogate 

#endif // SURROGATE_MODEL_DECLARATIONS 

