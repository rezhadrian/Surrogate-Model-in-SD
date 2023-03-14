/**
  * @file MasSpringDamper_IPCE.cpp 
  *
  * @brief 
  * Implementations of intrusive PCE model for Mass Spring Damper system 
  *
  * @anchor _MassSpringDamper_IPCE_cpp_ 
  * 
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "Surrogate_MassSpringDamper.hpp" 


namespace MassSpringDamper::Surrogate {

    IntrusivePCE::IntrusivePCE (

        const AnalyticalModel* SDModel, const R Omega, const Z Dim 

    ) : SDModel_( SDModel ), Omega_( Omega ), Dim_( Dim ) {}

} // Mass Spring Damper Intrusive PCE constructor 


namespace MassSpringDamper::Surrogate {

    void IntrusivePCE::SetIndices ( const Z MaxSum, const Z iMax ) {

        // Obtain sets of indices with sum of indices within the set <= MaxSum 
        auto Indices = BasisFunctions::MultiIndex<size_t> ( Dim(), MaxSum );

        // Remove sets with individual index > iMax 
        BasisFunctions::TotalTruncation<size_t> ( Indices, Dim(), iMax );

        Indices_ = Indices;

    }

} // Mass Spring Damper Intrusive PCE set indices 


namespace MassSpringDamper::Surrogate {

    void IntrusivePCE::Train (

        const VectorC& Load, 
        const VectorR& MassBasisCoeffs, 
        const VectorR& DamperBasisCoeffs, 
        const VectorR& SpringBasisCoeffs, 
        const VectorC& ForceBasisCoeffs 

    ) {

        auto nBasis = Indices_.size() / Dim_;
        auto nDOFs  = SDModel_ -> Dim ();

        auto DynamicStiffness = SDModel_ -> DynamicStiffness ( Omega_ );

        auto ExpHermiteTriples = BasisFunctions::ExpHermiteTriples<Z,R> (

            Indices_, 
            Dim(), 
            0

        );


        // ===================================================================
        // Calculate deterministic part of modified dynamic stiffness matrix 
        // ===================================================================

        Eigen::Map<MatrixXC, Eigen::RowMajor> dynamicStiffness ( 
            DynamicStiffness.data(), nDOFs, nDOFs  
        );

        Eigen::Map<MatrixXR> expHermiteTriples ( 
            ExpHermiteTriples.data(), nBasis, nBasis 
        );

        MatrixXC modDynamicStiffness = 
            MatrixXC::Zero ( nBasis * nDOFs, nBasis * nDOFs );

        for ( auto i = 0; i < nBasis; i++ ) {
        for ( auto j = 0; j < nBasis; j++ ) {

            modDynamicStiffness.block ( 

                i * nDOFs, 
                j * nDOFs, 

                nDOFs, nDOFs 

            ) += expHermiteTriples(i,j) * dynamicStiffness;

        }
        }


        // ===================================================================
        // Calculate random part of modified dynamic stiffness matrix 
        // ===================================================================

        auto nRandomBasis = MassBasisCoeffs.size() / nDOFs;

        for ( auto k = 0; k < nRandomBasis; k++ ) {

            VectorR RandomMasses ( 

                MassBasisCoeffs.begin() + k * nDOFs, 
                MassBasisCoeffs.begin() + k * nDOFs + nDOFs 

            );

            auto RandomMassMatrix = SDModel_ -> MassMatrix ( RandomMasses );

            VectorR RandomDampers ( 

                DamperBasisCoeffs.begin() + k * nDOFs, 
                DamperBasisCoeffs.begin() + k * nDOFs + nDOFs 

            );

            auto RandomDamping = SDModel_ -> DampingMatrix ( RandomDampers );

            VectorR RandomSprings ( 

                SpringBasisCoeffs.begin() + k * nDOFs, 
                SpringBasisCoeffs.begin() + k * nDOFs + nDOFs 

            );

            auto RandomStiffness = SDModel_ -> StiffnessMatrix ( RandomSprings );

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
            
            auto EHermiteTriples = BasisFunctions::ExpHermiteTriples<Z,R> (

                Indices_, 
                Dim(), 
                k

            );

            Eigen::Map<MatrixXR> eHermiteTriples ( 
                EHermiteTriples.data(), nBasis, nBasis 
            );

            for ( auto i = 0; i < nBasis; i++ ) {
            for ( auto j = 0; j < nBasis; j++ ) {

                modDynamicStiffness.block ( 

                    i * nDOFs, 
                    j * nDOFs, 

                    nDOFs, nDOFs 

                ) += eHermiteTriples(i,j) * randomDynStiffness;

            }
            }

        }


        // ===================================================================
        // Solve for coefficients of basis functions 
        // ===================================================================

        Coeffs_ = VectorC ( nDOFs * nBasis, 0.0 );

        VectorXC force = VectorXC::Zero( nDOFs * nBasis );

        Eigen::Map<const VectorXC> load ( 
            Load.data(), nDOFs 
        );

        Eigen::Map<VectorXC> coeffs ( 
            Coeffs_.data(), nDOFs * nBasis 
        );

        force.segment ( 0, nDOFs ) = load;

        Eigen::Map<const VectorXC> randomLoad ( 
            ForceBasisCoeffs.data(), ForceBasisCoeffs.size() 
        );

        force.segment ( 0, ForceBasisCoeffs.size() ) += randomLoad;

        coeffs = modDynamicStiffness.partialPivLu().solve(force);

    }

} // Mass Spring Damper Intrusive PCE train 

namespace MassSpringDamper::Surrogate {

    VectorC IntrusivePCE::ComputeResponse ( const VectorC& X ) const {

        auto nPoints = X.size() / Dim_;
        auto nBasis  = Indices_.size() / Dim_;
        auto nDOFs   = SDModel_ -> Dim ();

        // Compute basis functions with random inputs as arguments 
        auto Basis = BasisFunctions::HermitePolynomials<Z,R,C> (

                Indices_, X, Dim_ 

        );

        VectorC Response ( nDOFs * nPoints );


        // ===================================================================
        // Map vectors to eigen objects for linear algebra operations 
        // ===================================================================

        Eigen::Map<MatrixXC> basis ( 
            Basis.data(), nBasis , nPoints
        );

        Eigen::Map<const MatrixXC> coeffs (
            Coeffs_.data(), nDOFs, nBasis 
        );

        Eigen::Map<MatrixXC> response ( 
            Response.data(), nDOFs, nPoints 
        );


        // ===================================================================
        // Approximate response as linear combination of basis functions 
        // ===================================================================

        response = coeffs * basis;

        return Response;

    }

} // Mass Spring Damper Intrusive PCE compute response 




// namespace Surrogate {
//
//     void IntrusivePCE::PrintCoeffs () const {
//
//         typedef Eigen::Vector<size_t,Eigen::Dynamic> VectorZ;
//         Eigen::Map<const VectorXC> Coeffs ( Coeffs_.data(), Coeffs_.size() );
//         Eigen::Map<const VectorZ> Idx ( Indices_.data(), Indices_.size() );
//
//         std::cout << Coeffs << std::endl;
//
//     }
//
// } // Intrusive PCE PrintCoeffs 
//
//
// namespace Surrogate {
//
//     VectorC IntrusivePCE::Coeffs () const {
//
//         return Coeffs_;
//
//     }
//
// } // Intrusive PCE Coeffs 


