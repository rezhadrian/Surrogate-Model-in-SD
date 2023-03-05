/**
  * @file IntrusiveRPCE.cpp 
  *
  * @brief 
  * Implementations of intrusive RPCE model 
  *
  * @anchor _IntrusiveRPCE_cpp_ 
  * 
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "SurrogateModel.hpp" 
#include <ostream>


namespace Surrogate {

    IntrusiveRPCE::IntrusiveRPCE (

        const AnalyticalModel* SDModel, const R omega 

    ) : SDModel_( SDModel ), omega_( omega ), Dim_ ( SDModel_->Dim() ) {}

} // Intrusive RPCE constructor 


namespace Surrogate {

    void IntrusiveRPCE::SetNumIndices ( const Z iMax, const Z MaxSum ) {

        // Obtain sets of indices with max index <= iMax 
        auto NumIndices = BasisFunctions::MultiIndex<size_t> ( Dim(), iMax );

        // Remove sets with sum of indices > MaxSum 
        BasisFunctions::TotalTruncation<size_t> ( NumIndices, Dim(), MaxSum );

        NumIndices_ = NumIndices;

    }

} // Intrusive RPCE set numerator indices 


namespace Surrogate {

    void IntrusiveRPCE::SetDenIndices ( const Z iMax, const Z MaxSum ) {

        // Obtain sets of indices with max index <= iMax 
        auto DenIndices = BasisFunctions::MultiIndex<size_t> ( Dim(), iMax );

        // Remove sets with sum of indices > MaxSum 
        BasisFunctions::TotalTruncation<size_t> ( DenIndices, Dim(), MaxSum );

        DenIndices_ = VectorZ ( DenIndices.begin() + Dim (), DenIndices.end() );

    }

} // Intrusive RPCE set denominator indices 


namespace Surrogate {

    VectorC IntrusiveRPCE::ComputeResponse ( const VectorC& X ) const {

        auto nPoints = X.size() / Dim();


        // ===================================================================
        // Compute numerator to approximate response 
        // ===================================================================

        auto NumBasis = BasisFunctions::HermitePolynomials<Z,R,C> (
                NumIndices_, X, Dim() 
        );

        VectorC Num ( nPoints * Dim() );

        Eigen::Map<MatrixXC> numBasis (
            NumBasis.data(), 
            NumIndices_.size() / Dim(), 
            nPoints 
        );

        Eigen::Map<const MatrixXC> numCoeffs (
            NumCoeffs_.data(), Dim(), NumIndices_.size() / Dim() 
        );

        Eigen::Map<MatrixXC> num ( 
            Num.data(), Dim(), nPoints 
        );

        num = numCoeffs * numBasis;


        // ===================================================================
        // Compute denominator to approximate response 
        // ===================================================================

        auto DenBasis = BasisFunctions::HermitePolynomials<Z,R,C> (
                DenIndices_, X, Dim()
        );

        VectorC Den ( nPoints * Dim() );

        Eigen::Map<MatrixXC> denBasis (
            DenBasis.data(), 
            DenIndices_.size() / Dim(), 
            nPoints 
        );

        Eigen::Map<const MatrixXC> denCoeffs (
            DenCoeffs_.data(), Dim(), DenIndices_.size() / Dim() 
        );

        Eigen::Map<MatrixXC> den ( 
            Den.data(), Dim(), nPoints 
        );

        den = denCoeffs * denBasis;


        // ===================================================================
        // Approximate response as numerator / ( 1 + denominator ) 
        // ===================================================================

        VectorC Response ( nPoints * Dim(), 0.0 );

        std::transform (

            Num.begin(), Num.end(), 
            Den.begin(), 

            Response.begin(), 

            [](const auto numerator, const auto denominator) {

                return numerator / ( C(1.0) + denominator );
            }
            
        );

        return Response;

    }

} // Intrusive RPCE compute response 


namespace Surrogate {

    void IntrusiveRPCE::Train (
        
        const VectorC& Load 

    ) {


        // ===================================================================
        // Obtain minimum required number of basis functions in PCE model 
        // ===================================================================

        size_t iMax = 1;
        size_t MaxSum = 1;

        VectorZ Indices;

        while ( true ) {

            auto indices = BasisFunctions::MultiIndex<Z> ( Dim(), iMax );

            BasisFunctions::TotalTruncation<Z> ( indices, Dim(), MaxSum );

            if ( indices.size() > NumIndices_.size() + DenIndices_.size() ) {
                Indices = indices;
                break;
            }

            iMax++;
            MaxSum++;

        }


        // ===================================================================
        // Obtain coefficients of basis functions in PCE model 
        // ===================================================================

        IntrusivePCE PCEModel ( SDModel_, omega_ );

        PCEModel.SetIndices ( iMax, MaxSum );

        PCEModel.Train ( Load );

        auto PCECoeffs = PCEModel.Coeffs();


        // ===================================================================
        // Compute coefficients of basis functions in denominator 
        // ===================================================================

        auto nPCEBasis =     Indices.size() / Dim();
        auto nDenBasis = DenIndices_.size() / Dim();
        auto nNumBasis = NumIndices_.size() / Dim();

        MatrixXC modDenStiffness = MatrixXC::Zero ( 
            nDenBasis * Dim(), 
            nDenBasis * Dim() 
        );

        VectorXC modDenForce = VectorXC::Zero ( 
            nDenBasis * Dim() 
        );

        DenCoeffs_ = VectorC ( nDenBasis * Dim(), 0.0 );

        Eigen::Map<MatrixXC> pceCoeffs ( 
            PCECoeffs.data(), Dim(), PCECoeffs.size() / Dim() 
        );

        Eigen::Map<VectorXC> denCoeffs ( 
            DenCoeffs_.data(), nDenBasis * Dim() 
        );

        for ( auto j = 1; j < nDenBasis+1; j++ ) {

            auto ExpHermiteTriples = BasisFunctions::ExpHermiteTriples<Z,R> (
                    Indices, 
                    Dim(), 
                    j
            );

            Eigen::Map<MatrixXR> expHermiteTriples ( 
                ExpHermiteTriples.data(), nPCEBasis, nPCEBasis 
            );

            auto YHT = pceCoeffs * expHermiteTriples;

            for ( auto l = 1; l < nDenBasis+1; l++ ) {

                modDenStiffness.block(

                    (l-1) * Dim(),
                    (j-1) * Dim(),

                    Dim(),
                    Dim()

                ) = YHT.col(nNumBasis+l).asDiagonal();

                modDenForce.segment ( 

                    (l-1) * Dim(), 

                    Dim() 

                ) = -pceCoeffs.col(nNumBasis+l);

            }

        }

        denCoeffs = modDenStiffness.fullPivLu().solve(modDenForce);


        // ===================================================================
        // Compute coefficients of basis functions in numerator  
        // ===================================================================

        MatrixXC modNumStiffness = MatrixXC::Zero ( 
            nNumBasis * Dim(), 
            nDenBasis * Dim() 
        );

        VectorXC modNumForce = VectorXC::Zero ( 
            nNumBasis * Dim() 
        );

        NumCoeffs_ = VectorC ( nNumBasis * Dim(), 0.0 );

        Eigen::Map<VectorXC> numCoeffs ( 
            NumCoeffs_.data(), nNumBasis * Dim() 
        );

        for ( auto j = 1; j < nDenBasis+1; j++ ) {

            auto ExpHermiteTriples = BasisFunctions::ExpHermiteTriples<Z,R> (
                Indices, 
                Dim(), 
                j
            );

            Eigen::Map<MatrixXR> expHermiteTriples ( 
                ExpHermiteTriples.data(), nPCEBasis, nPCEBasis 
            );

            auto YHT = pceCoeffs * expHermiteTriples;

            for ( auto l = 0; l < nNumBasis; l++ ) {

                modNumStiffness.block(

                    l * Dim(),
                    (j-1) * Dim(),
                    Dim(),
                    Dim()

                ) = YHT.col(l).asDiagonal();

                modNumForce.segment ( 

                    l * Dim(), 

                    Dim() 

                ) = pceCoeffs.col(l);

            }

        }

        numCoeffs = modNumStiffness * denCoeffs + modNumForce;

    }

} // Intrusive RPCE train 

