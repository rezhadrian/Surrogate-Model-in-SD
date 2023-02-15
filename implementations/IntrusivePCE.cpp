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

        // X contains number of points x dimension of each point 
        auto nPoints = X.size() / Dim ();

        // Compute basis functions with X as arguments 
        auto Basis = 
            BasisFunctions::HermitePolynomials<Z,R,C> (

                Indices_, X, Dim ()

        );

        // Map vectors to eigen object for linear algebra operations 

        Eigen::Map<MatrixXC> Psi ( 
            Basis.data(), Indices_.size() / Dim (), nPoints
        );

        Eigen::Map<const MatrixXC> Coeffs (
            Coeffs_.data(), Dim (), Indices_.size() / Dim ()
        );

        VectorC result ( nPoints * Dim () );

        Eigen::Map<MatrixXC> Result ( result.data(), Dim (), nPoints );

        // std::cout << Psi << std::endl; 

        Result = Coeffs * Psi;

        return result;

    }

} // Intrusive PCE compute response 


namespace Surrogate {

    void IntrusivePCE::Train (

        const VectorC& Load

    ) {

        auto P = Indices_.size() / Dim ();

        MatrixXC H = MatrixXC::Zero ( P * Dim(), P * Dim() );

        auto KD = SDModel_ -> DynamicStiffness ( omega_ );
        // VectorC KD ( Dim () );

        auto HermiteTriples = BasisFunctions::ExpHermiteTriples<Z,R> (

            Indices_, 
            Dim(), 
            0

        );

        Eigen::Map<MatrixXC, Eigen::RowMajor> KDD ( KD.data(), Dim(), Dim() );
        Eigen::Map<MatrixXR> HT ( HermiteTriples.data(), P, P );

        // std::cout << HT << std::endl;

        for ( auto i = 0; i < P; i++ ) {
        for ( auto j = 0; j < P; j++ ) {

            H.block ( i*Dim(), j*Dim(), Dim(), Dim() ) +=
                HT(i,j) * KDD;

        }
        }

        for ( auto k = 0; k < Dim(); k++ ) {

            VectorR Kkv ( Dim(), 0.0 );
            Kkv[k] = 1.0;

            auto Kks = SDModel_ -> StiffnessMatrix ( Kkv );
            // VectorR Kks ( Dim () );

            auto HermiteTriples = BasisFunctions::ExpHermiteTriples<Z,R> (

                Indices_, 
                Dim(), 
                k+1

            );

            // for ( auto ii = 0; ii < Indices_.size(); ii++ ) {
            //
            //     std::cout << Indices_[ii] << std::endl;
            //    
            // }

            

            Eigen::Map<MatrixXR, Eigen::RowMajor> KDD ( Kks.data(), Dim(), Dim() );
            Eigen::Map<MatrixXR> HT ( HermiteTriples.data(), P, P );

            // std::cout << HT << std::endl;
            // std::cout << KDD << std::endl;

            for ( auto i = 0; i < P; i++ ) {
            for ( auto j = 0; j < P; j++ ) {

                H.block ( i*Dim(), j*Dim(), Dim(), Dim() ) +=
                    HT(i,j) * KDD;

            }
            }
        }

        // std::cout << H << std::endl;

        Coeffs_ = VectorC ( P * Dim (), 0.0 );

        Eigen::Map<const VectorXC> Force ( Load.data(), Dim() );

        VectorXC F = VectorXC::Zero( P * Dim () );
        F.segment ( 0, Dim() ) = Force;

        Eigen::Map<VectorXC> Coeffs ( Coeffs_.data(), P * Dim () );

        Coeffs = H.ldlt().solve(F);

    }

} // Intrusive PCE train 


namespace Surrogate {

    void IntrusivePCE::PrintCoeffs () const {

        typedef Eigen::Vector<size_t,Eigen::Dynamic> VectorZ;
        Eigen::Map<const VectorXC> Coeffs ( Coeffs_.data(), Coeffs_.size() );
        Eigen::Map<const VectorZ> Idx ( Indices_.data(), Indices_.size() );

        std::cout << Coeffs << std::endl;

    }

}


