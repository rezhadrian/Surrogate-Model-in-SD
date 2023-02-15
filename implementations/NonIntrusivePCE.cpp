/**
  * @file NonIntrusivePCE.cpp 
  *
  * @brief 
  * Implementations of non-intrusive PCE model 
  *
  * @anchor _NonIntrusivePCE_cpp_ 
  * 
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "SurrogateModel.hpp" 


namespace Surrogate {

    NonIntrusivePCE::NonIntrusivePCE ( 

        const AnalyticalModel* SDModel, const R omega 

    ) : SDModel_( SDModel ), omega_( omega ), Dim_(SDModel_->Dim()) {}

} // Non-intrusive PCE constructor 


namespace Surrogate {

    void NonIntrusivePCE::SetIndices ( Z iMax, Z MaxSum ) {

        // Obtain sets of indices with max index <= iMax 
        auto Indices = BasisFunctions::MultiIndex<size_t> ( Dim(), iMax );

        // Remove sets with sum of indices > MaxSum 
        BasisFunctions::TotalTruncation<size_t> ( Indices, Dim(), MaxSum );

        Indices_ = Indices;

    }

} // Non-intrusive PCE set indices 


namespace Surrogate {

    VectorC NonIntrusivePCE::ComputeResponse ( const VectorC& X ) const {

        // X contains number of points x dimension of each point 
        auto nPoints = X.size() / Dim ();

        // Compute basis functions with X as arguments 
        auto Basis = 
            BasisFunctions::HermitePolynomials<Z,R,C> (

                Indices_, X, Dim ()

        );

        // Map vectors to eigen object for linear algebra operations 

        Eigen::Map<MatrixXC, Eigen::RowMajor> Psi ( 
            Basis.data(), nPoints, Indices_.size() / Dim ()
        );

        Eigen::Map<const MatrixXC> Coeffs (
            Coeffs_.data(), Indices_.size() / Dim (), Dim ()
        );

        VectorC result ( nPoints * Dim () );

        Eigen::Map<MatrixXC> Result ( result.data(), nPoints, Dim () );

        Result = Psi * Coeffs;

        return result;

    }

} // Non-intrusive PCE ComputeResponse 


namespace Surrogate {

    void NonIntrusivePCE::Train ( 

        const VectorC& Load, const VectorR& TrainSet 

    ) {

        auto nPoints = TrainSet.size() / Dim ();

        VectorC Response ( TrainSet.size() );

        // Compute response for each point 
        for ( auto i = 0; i < nPoints; i++ ) {

            auto start = TrainSet.begin() + i * Dim ();
            auto end   = TrainSet.begin() + i * Dim () + Dim ();

            auto Disp = SDModel_ -> ComputeResponse ( 
                Load, omega_, start, end 
            );

            auto source_start = Disp.begin(); 
            auto source_end   = Disp.end(); 

            auto dest_start   = Response.begin() + i * Dim ();
            auto dest_end     = Response.begin() + i * Dim () + Dim ();

            Response.erase ( dest_start, dest_end );
            Response.insert ( 
                dest_start, 
                std::make_move_iterator ( source_start ), 
                std::make_move_iterator ( source_end   )
            );

        }

        VectorC TSet ( TrainSet.size(), 0.0 );

        std::transform ( 

            TrainSet.begin(), TrainSet.end(), 
            TSet.begin(), 

            [](const auto m){return C(m);}

        );

        auto Basis = 
            BasisFunctions::HermitePolynomials<Z,R,C> (

                Indices_, TSet, Dim ()

        );

        Coeffs_ = VectorC ( Indices_.size() , 0.0 );

        Eigen::Map<MatrixXC, Eigen::RowMajor> Psi ( 
            Basis.data(), nPoints, Indices_.size() / Dim ()
        );

        Eigen::Map<MatrixXC> Coeffs (
            Coeffs_.data(), Indices_.size() / Dim (), Dim ()
        );

        Eigen::Map<MatrixXC, Eigen::RowMajor> Result (
            Response.data(), nPoints, Dim ()
        );

        Coeffs = (Psi.transpose() * Psi).ldlt().solve ( 
                Psi.transpose() * Result );

    }

} // Non-intrusive PCE Train 


namespace Surrogate {

    void NonIntrusivePCE::PrintCoeffs () const {

        typedef Eigen::Vector<size_t,Eigen::Dynamic> VectorZ;
        Eigen::Map<const VectorXC> Coeffs ( Coeffs_.data(), Coeffs_.size() );
        Eigen::Map<const VectorZ> Idx ( Indices_.data(), Indices_.size() );

        std::cout << Coeffs << std::endl;
        std::cout << Idx << std::endl;

    }

} // Non-intrusive PCE PrintCoeffs 

