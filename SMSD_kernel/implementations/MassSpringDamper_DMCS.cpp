/**
  * @file MasSpringDamper_DMCS.cpp 
  *
  * @brief 
  * Implementations direct MCS class for Mass Spring Damper system 
  *
  * @anchor _MassSpringDamper_DMCS_cpp_ 
  * 
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "Surrogate_MassSpringDamper.hpp" 


namespace MassSpringDamper::Surrogate {

    DirectMCS::DirectMCS (

        const VectorR& Masses, 
        const VectorR& Dampers, 
        const VectorR& Springs, 
        const R Omega, 
        const Z Dim 

    ) : 
        Masses_ ( Masses ), 
        Dampers_ ( Dampers ), 
        Springs_ ( Springs ), 
        Omega_( Omega ), Dim_( Dim ) {}

} // Mass Spring Damper direct MCS constructor 


namespace MassSpringDamper::Surrogate {

    void DirectMCS::SetIndices ( const Z MaxSum, const Z iMax ) {

        // Obtain sets of indices with sum of indices within the set <= MaxSum 
        auto Indices = BasisFunctions::MultiIndex<size_t> ( Dim_, MaxSum );

        // Remove sets with individual index > iMax 
        BasisFunctions::TotalTruncation<size_t> ( Indices, Dim_, iMax );

        Indices_ = Indices;

    }

} // Mass Spring Damper direct MCS set indices 


namespace MassSpringDamper::Surrogate {

    VectorC DirectMCS::ComputeResponse (

        const VectorC& X, 
        const VectorC& Load, 
        const VectorR& MassBasisCoeffs, 
        const VectorR& DamperBasisCoeffs, 
        const VectorR& SpringBasisCoeffs, 
        const VectorC& ForceBasisCoeffs 

    ) const {

        AnalyticalModel SDModel_ ( Masses_, Dampers_, Springs_ );

        auto nPoints = X.size() / Dim_; 
        auto nBasis  = Indices_.size() / Dim_; 
        auto nDOFs   = SDModel_.Dim ();

        auto Basis = BasisFunctions::HermitePolynomials<Z,R,C> (

            Indices_, X, Dim_ 

        );

        Eigen::Map <MatrixXC> basis ( 
            Basis.data(), nBasis, nPoints 
        );

        VectorR RandomMasses ( nDOFs * nPoints ); 

        Eigen::Map <MatrixXR> randomMasses (
            RandomMasses.data(), nDOFs, nPoints 
        );

        Eigen::Map <const MatrixXR> massBasisCoeffs (
            MassBasisCoeffs.data(), nDOFs, nBasis 
        );

        randomMasses = massBasisCoeffs * basis.real(); 

        
        VectorR RandomDampers ( nDOFs * nPoints ); 

        Eigen::Map <MatrixXR> randomDampers (
            RandomDampers.data(), nDOFs, nPoints 
        );

        Eigen::Map <const MatrixXR> damperBasisCoeffs (
            DamperBasisCoeffs.data(), nDOFs, nBasis 
        );

        randomDampers = damperBasisCoeffs * basis.real(); 


        VectorR RandomSprings ( nDOFs * nPoints ); 

        Eigen::Map <MatrixXR> randomSprings (
            RandomSprings.data(), nDOFs, nPoints 
        );

        Eigen::Map <const MatrixXR> springBasisCoeffs (
            SpringBasisCoeffs.data(), nDOFs, nBasis 
        );

        randomSprings = springBasisCoeffs * basis.real(); 


        VectorC RandomForces ( nDOFs * nPoints );

        Eigen::Map <MatrixXC> randomForces (
            RandomForces.data(), nDOFs, nPoints 
        );

        Eigen::Map <const MatrixXC> forceBasisCoeffs (
            ForceBasisCoeffs.data(), nDOFs, nBasis 
        );

        randomForces = forceBasisCoeffs * basis; 


        VectorC Result ( nDOFs * nPoints );

        Eigen::Map <MatrixXC> result (
            Result.data(), nDOFs, nPoints 
        );

        for ( auto i = 0; i < nPoints; i++ ) {

            VectorR Masses ( nDOFs );

            std::transform (

                RandomMasses.begin() + i * nDOFs, 
                RandomMasses.begin() + i * nDOFs + nDOFs, 
                Masses_.begin(), 

                Masses.begin(), 

                [](const auto m, const auto n) { return m+n; }

            );

            auto RandomMassMatrix = SDModel_.MassMatrix ( Masses );


            VectorR Dampers ( nDOFs );

            std::transform (

                RandomDampers.begin() + i * nDOFs, 
                RandomDampers.begin() + i * nDOFs + nDOFs, 
                Dampers_.begin(), 

                Dampers.begin(), 

                [](const auto m, const auto n) { return m+n; }

            );

            auto RandomDamping = SDModel_.DampingMatrix ( Dampers );

            VectorR Springs ( nDOFs );

            std::transform (

                RandomSprings.begin() + i * nDOFs, 
                RandomSprings.begin() + i * nDOFs + nDOFs, 
                Springs_.begin(), 

                Springs.begin(), 

                [](const auto m, const auto n) { return m+n; }

            );

            auto RandomStiffness = SDModel_.StiffnessMatrix ( Springs );


            VectorC Forces ( nDOFs );

            std::transform (

                RandomForces.begin() + i * nDOFs, 
                RandomForces.begin() + i * nDOFs + nDOFs, 
                Load.begin(), 

                Forces.begin(), 

                [](const auto m, const auto n) { return m+n; }

            );

            VectorC RandomDynStiffness ( nDOFs * nDOFs, 0.0 );

            Eigen::Map <MatrixXR> randomMassMatrix ( 
                RandomMassMatrix.data(), nDOFs, nDOFs 
            );

            Eigen::Map <MatrixXR> randomDamping ( 
                RandomDamping.data(), nDOFs, nDOFs 
            );

            Eigen::Map <MatrixXR> randomStiffness ( 
                RandomStiffness.data(), nDOFs, nDOFs 
            );

            Eigen::Map <MatrixXC> randomDynStiffness ( 
                RandomDynStiffness.data(), nDOFs, nDOFs 
            );

            randomDynStiffness.real() = 
                randomStiffness - Omega_ * Omega_ * randomMassMatrix;

            randomDynStiffness.imag() = 
                Omega_ * randomDamping;

            Eigen::Map <VectorXC> forces (
                Forces.data(), nDOFs 
            );

            result.col(i) = randomDynStiffness.partialPivLu().solve(forces);

        }

        return Result;

    }


}


