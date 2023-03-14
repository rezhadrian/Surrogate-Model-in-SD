/**
  * @file Surrogate_MassSpringDamper.hpp 
  * 
  * @brief 
  * Declarations of Surrogate Models for Mass Spring Damper System 
  * 
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef SURROGATE_MASS_SPRING_DAMPER_DECLARATIONS 
#define SURROGATE_MASS_SPRING_DAMPER_DECLARATIONS 

#include "BasisFunctions.hpp" 
#include "MonteCarlo.hpp" 
#include "AnalyticalModel.hpp" 
#include "LibrariesLoader_SM.hpp" 


namespace MassSpringDamper::Surrogate {


    typedef size_t          Z;
    typedef double          R;
    typedef std::complex<R> C; 

    typedef std::vector<R> VectorR; 
    typedef std::vector<C> VectorC; 
    typedef std::vector<Z> VectorZ; 

    typedef Analytical::MassSpringDamper<Z,R,C> AnalyticalModel;

    typedef Eigen::Matrix<R, Eigen::Dynamic, Eigen::Dynamic> MatrixXR;
    typedef Eigen::Matrix<C, Eigen::Dynamic, Eigen::Dynamic> MatrixXC;

    typedef Eigen::Vector<C, Eigen::Dynamic> VectorXC;

    /**
      * @class DirectMCS 
      * 
      * @brief 
      * Class to perform direct MCS for Mass Spring Damper system. 
      */
    class DirectMCS {

        VectorR Masses_; 
        VectorR Dampers_; 
        VectorR Springs_; 

        VectorZ Indices_; 

        R Omega_; 
        Z Dim_; 

        public: 

        /**
          * @brief 
          * Create PCE model for given SD model and angular velocity 
          * 
          * @param SDModel address of Mass Spring Damper model 
          * @param omega   angular velocity of harmonic load 
          */
        DirectMCS ( 

            const VectorR& Masses, 
            const VectorR& Dampers, 
            const VectorR& Springs, 
            const R Omega, 
            const Z Dim 

        );

        /**
          * @brief 
          * Set indices for PCE basis functions 
          * 
          * @param iMax   largest allowable sum of indices in a set 
          * @param MaxSum largest allowable individual index 
          */
        void SetIndices ( const Z iMax, const Z MaxSum );

        /**
          * @brief 
          * Compute coefficients of basis functions 
          * 
          * @param Load harmonic load vector 
          */
        VectorC ComputeResponse ( 

            const VectorC& X, 
            const VectorC& Load, 
            const VectorR& MassBasisCoeffs, 
            const VectorR& DamperBasisCoeffs, 
            const VectorR& SpringBasisCoeffs, 
            const VectorC& ForceBasisCoeffs 

        ) const; 

    }; // DirectMCS 

    /**
      * @class IntrusivePCE 
      * 
      * @brief 
      * Intrusive surrogate model for Mass Spring Damper system. 
      * Use weak formulation based on polynomial chaos expansion. 
      * Implemented in @ref _MassSpringDamper_IPCE_cpp_ 
      */
    class IntrusivePCE {

        const AnalyticalModel* SDModel_;

        VectorZ Indices_; 
        VectorC Coeffs_; 

        R Omega_; 
        Z Dim_; 

        public: 

        Z Dim () const { return Dim_; } 

        /**
          * @brief 
          * Create PCE model for given SD model and angular velocity 
          * 
          * @param SDModel address of Mass Spring Damper model 
          * @param omega   angular velocity of harmonic load 
          */
        IntrusivePCE ( 

            const AnalyticalModel* SDModel, 
            const R Omega, 
            const Z Dim 

        );

        /**
          * @brief 
          * Set indices for PCE basis functions 
          * 
          * @param iMax   largest allowable sum of indices in a set 
          * @param MaxSum largest allowable individual index 
          */
        void SetIndices ( const Z iMax, const Z MaxSum );

        /**
          * @brief 
          * Compute coefficients of basis functions 
          * 
          * @param Load harmonic load vector 
          */
        void Train ( 

            const VectorC& Load, 
            const VectorR& MassBasisCoeffs, 
            const VectorR& DamperBasisCoeffs, 
            const VectorR& SpringBasisCoeffs, 
            const VectorC& ForceBasisCoeffs 

        ); 

        /**
          * @brief 
          * Approximate response of analytical model for a given random inputs 
          * 
          * @param Load deterministic load 
          * @param MassBasisCoeffs   coefficients of PCEs added to masses 
          * @param DamperBasisCoeffs coefficients of PCEs added to dampers 
          * @param SpringBasisCoeffs coefficients of PCEs added to springs 
          * @param ForceBasisCoeffs  coefficients of PCEs added to forces 
          * 
          * @return approximate displacement vector 
          */
        VectorC ComputeResponse ( const VectorC& Load ) const;

    }; // IntrusivePCE 


} // MassSpringDamper::Surrogate 


#endif // SURROGATE_MASS_SPRING_DAMPER_DECLARATIONS 

