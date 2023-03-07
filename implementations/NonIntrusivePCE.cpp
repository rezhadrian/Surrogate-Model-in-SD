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
#include <iostream>


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

        Eigen::Map<MatrixXC> Psi ( 
            Basis.data(), Indices_.size() / Dim (), nPoints
        );

        Eigen::Map<const MatrixXC> Coeffs (
            Coeffs_.data(), Indices_.size() / Dim (), Dim ()
        );

        VectorC result ( nPoints * Dim () );

        // Eigen::Map<MatrixXC> Result ( result.data(), nPoints, Dim () );
        Eigen::Map<MatrixXC> Result ( 
            result.data(),  Dim(), nPoints 
        );

        Result = Coeffs.transpose() * Psi;

        // std::cout << Result << std::endl;

        // for ( auto i = 0; i < result.size(); i++ ) {
        //     std::cout << result[i] << std::endl;
        // }

        return result;

    }

} // Non-intrusive PCE ComputeResponse 


namespace Surrogate {

    void NonIntrusivePCE::Train ( 

        const VectorC& Load, const VectorR& TrainSet 

    ) {

        auto nPoints = TrainSet.size() / Dim ();

        // for ( auto i = 0; i < TrainSet.size(); i++ ) {
        //     std::cout << TrainSet[i] << std::endl;
        // }

        VectorC Response ( TrainSet.size() );

        // Compute response for each point 
        for ( auto i = 0; i < nPoints; i++ ) {

            auto start = TrainSet.begin() + i * Dim ();
            auto end   = TrainSet.begin() + i * Dim () + Dim ();

            auto Disp = SDModel_ -> ComputeResponse ( 
                Load, omega_, start, end 
            );

            // for ( auto ii = 0; ii < Disp.size(); ii++ ) {
            //
            //     std::cout << Disp[ii] << std::endl;
            //
            // }

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

        Eigen::Map<MatrixXC> Psi ( 
            Basis.data(), Indices_.size() / Dim (), nPoints 
        );

        // std::cout << Psi << std::endl << std::endl;

        Eigen::Map<MatrixXC> Coeffs (
            Coeffs_.data(), Indices_.size() / Dim (), Dim ()
        );

        Eigen::Map<MatrixXC> Result (
            Response.data(), Dim (), nPoints 
        );

        // std::cout << Result << std::endl << std::endl;;

        // std::cout << Result << std::endl;

        Coeffs = (Psi * Psi.transpose()).partialPivLu().solve ( 
                Psi * Result.transpose() );

        // std::cout << Coeffs << std::endl << std::endl;

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

