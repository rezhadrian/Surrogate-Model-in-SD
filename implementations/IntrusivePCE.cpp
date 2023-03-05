/**
  * @file IntrusivePCE.cpp 
  *
  * @brief 
  * Implementations of intrusive PCE model 
  *
  * @anchor _IntrusivePCE_cpp_ 
  * 
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "SurrogateModel.hpp" 


namespace Surrogate {

    IntrusivePCE::IntrusivePCE (

        const AnalyticalModel* SDModel, const R omega 

    ) : SDModel_( SDModel ), omega_( omega ), Dim_( SDModel_->Dim() ) {}

} // Intrusive PCE constructor 


namespace Surrogate {

    void IntrusivePCE::SetIndices ( Z iMax, Z MaxSum ) {

        // Obtain sets of indices with max index <= iMax 
        auto Indices = BasisFunctions::MultiIndex<size_t> ( Dim(), iMax );

        // Remove sets with sum of indices > MaxSum 
        BasisFunctions::TotalTruncation<size_t> ( Indices, Dim(), MaxSum );

        Indices_ = Indices;

    }

} // Intrusive PCE set indices 


namespace Surrogate {

    VectorC IntrusivePCE::ComputeResponse ( const VectorC& X ) const {

        auto nPoints = X.size() / Dim();

        // Compute basis functions with random inputs as arguments 
        auto Basis = BasisFunctions::HermitePolynomials<Z,R,C> (

                Indices_, X, Dim()

        );

        VectorC Response ( nPoints * Dim () );


        // ===================================================================
        // Map vectors to eigen objects for linear algebra operations 
        // ===================================================================

        Eigen::Map<MatrixXC> basis ( 
            Basis.data(), Indices_.size() / Dim (), nPoints
        );

        Eigen::Map<const MatrixXC> coeffs (
            Coeffs_.data(), Dim (), Indices_.size() / Dim ()
        );

        Eigen::Map<MatrixXC> response ( 
            Response.data(), Dim (), nPoints 
        );


        // ===================================================================
        // Approximate response as linear combination of basis functions 
        // ===================================================================

        response = coeffs * basis;

        return Response;

    }

} // Intrusive PCE compute response 


namespace Surrogate {

    void IntrusivePCE::Train (

        const VectorC& Load  

    ) {

        auto nBasis = Indices_.size() / Dim();

        auto DynamicStiffness = SDModel_ -> DynamicStiffness ( omega_ );

        auto ExpHermiteTriples = BasisFunctions::ExpHermiteTriples<Z,R> (

            Indices_, 
            Dim(), 
            0

        );


        // ===================================================================
        // Calculate deterministic part of modified dynamic stiffness matrix 
        // ===================================================================

        Eigen::Map<MatrixXC, Eigen::RowMajor> dynamicStiffness ( 
            DynamicStiffness.data(), Dim(), Dim() 
        );

        Eigen::Map<MatrixXR> expHermiteTriples ( 
            ExpHermiteTriples.data(), nBasis, nBasis 
        );

        MatrixXC modDynamicStiffness = 
            MatrixXC::Zero ( nBasis * Dim(), nBasis * Dim() );

        for ( auto i = 0; i < nBasis; i++ ) {
        for ( auto j = 0; j < nBasis; j++ ) {

            modDynamicStiffness.block ( 

                i*Dim(), 
                j*Dim(), 

                Dim(), Dim() 

            ) += expHermiteTriples(i,j) * dynamicStiffness;

        }
        }


        // ===================================================================
        // Calculate random part of modified dynamic stiffness matrix 
        // ===================================================================

        for ( auto k = 0; k < Dim(); k++ ) {

            VectorR RandomSprings ( Dim(), 0.0 );
            RandomSprings[k] = 1.0;

            auto RandomStiffness = SDModel_->StiffnessMatrix ( RandomSprings );

            auto EHermiteTriples = BasisFunctions::ExpHermiteTriples<Z,R> (

                Indices_, 
                Dim(), 
                k+1

            );

            Eigen::Map<MatrixXR, Eigen::RowMajor> randomStiffness ( 
                RandomStiffness.data(), Dim(), Dim() 
            );

            Eigen::Map<MatrixXR> eHermiteTriples ( 
                EHermiteTriples.data(), nBasis, nBasis 
            );

            for ( auto i = 0; i < nBasis; i++ ) {
            for ( auto j = 0; j < nBasis; j++ ) {

                modDynamicStiffness.block ( 

                    i*Dim(), 
                    j*Dim(), 

                    Dim(), Dim() 

                ) += eHermiteTriples(i,j) * randomStiffness;

            }
            }

        }


        // ===================================================================
        // Solve for coefficients of basis functions 
        // ===================================================================

        Coeffs_ = VectorC ( nBasis * Dim(), 0.0 );

        VectorXC force = VectorXC::Zero( nBasis * Dim() );

        Eigen::Map<const VectorXC> load ( 
            Load.data(), Dim() 
        );

        Eigen::Map<VectorXC> coeffs ( 
            Coeffs_.data(), nBasis * Dim () 
        );

        force.segment ( 0, Dim() ) = load;

        coeffs = modDynamicStiffness.partialPivLu().solve(force);

    }

} // Intrusive PCE train 


namespace Surrogate {

    void IntrusivePCE::PrintCoeffs () const {

        typedef Eigen::Vector<size_t,Eigen::Dynamic> VectorZ;
        Eigen::Map<const VectorXC> Coeffs ( Coeffs_.data(), Coeffs_.size() );
        Eigen::Map<const VectorZ> Idx ( Indices_.data(), Indices_.size() );

        std::cout << Coeffs << std::endl;

    }

} // Intrusive PCE PrintCoeffs 


namespace Surrogate {

    VectorC IntrusivePCE::Coeffs () const {

        return Coeffs_;

    }

} // Intrusive PCE Coeffs 


