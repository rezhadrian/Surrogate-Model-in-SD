/**
  * @file NonIntrusiveRPCE.cpp 
  *
  * @brief 
  * Implementations of non-intrusive RPCE model 
  *
  * @anchor _NonIntrusiveRPCE_cpp_ 
  * 
  * @author 
  * Rezha Adrian Tanuharja @n 
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#include "SurrogateModel.hpp" 


namespace Surrogate {

    NonIntrusiveRPCE::NonIntrusiveRPCE ( 

        const AnalyticalModel* SDModel, const R omega 

    ) : SDModel_( SDModel ), omega_( omega ), Dim_(SDModel_->Dim()) {}

} // Non-intrusive RPCE constructor 


namespace Surrogate {

    void NonIntrusiveRPCE::SetNumIndices ( Z iMax, Z MaxSum ) {

        // Obtain sets of indices with max index <= iMax 
        auto Indices = BasisFunctions::MultiIndex<size_t> ( Dim(), iMax );

        // Remove sets with sum of indices > MaxSum 
        BasisFunctions::TotalTruncation<size_t> ( Indices, Dim(), MaxSum );

        NumIndices_ = Indices;

    }

} // Non-intrusive RPCE SetNumIndices 


namespace Surrogate {

    void NonIntrusiveRPCE::SetDenIndices ( Z iMax, Z MaxSum ) {

        // Obtain sets of indices with max index <= iMax 
        auto Indices = BasisFunctions::MultiIndex<size_t> ( Dim(), iMax );

        // Remove sets with sum of indices > MaxSum 
        BasisFunctions::TotalTruncation<size_t> ( Indices, Dim(), MaxSum );

        DenIndices_ = Indices;

    }

} // Non-intrusive RPCE SetDenIndices 


namespace Surrogate {

    VectorC NonIntrusiveRPCE::ComputeResponse ( const VectorC& X ) const {


        // X contains number of points x dimension of each point 
        auto nPoints = X.size() / Dim ();

        // Compute numerator basis functions with X as arguments 
        auto NumBasis = 
            BasisFunctions::HermitePolynomials<Z,R,C> (

                NumIndices_, X, Dim ()

        );

        // Compute numerator basis functions with X as arguments 
        auto DenBasis = 
            BasisFunctions::HermitePolynomials<Z,R,C> (

                DenIndices_, X, Dim ()

        );

        Eigen::Map<MatrixXC> NumPsi ( 
            NumBasis.data(), NumIndices_.size() / Dim (), nPoints
        );

        Eigen::Map<MatrixXC> DenPsi ( 
            DenBasis.data(), DenIndices_.size() / Dim (), nPoints
        );

        Eigen::Map<const MatrixXC> NumCoeffs (
            NumCoeffs_.data(), NumIndices_.size() / Dim (), Dim ()
        );

        Eigen::Map<const MatrixXC> DenCoeffs (
            DenCoeffs_.data(), DenIndices_.size() / Dim (), Dim ()
        );

        VectorC result ( nPoints * Dim () );

        Eigen::Map<MatrixXC> Result ( result.data(), Dim (), nPoints );

        Result = ( NumCoeffs.transpose() * NumPsi ).array() / 
                 ( DenCoeffs.transpose() * DenPsi ).array();

        // std::cout << Result << std::endl;

        // for ( auto i = 0; i < result.size(); i++ ) {
        //     std::cout << result[i] << std::endl;
        // }

        return result;

    }

} // Non-intrusive RPCE ComputeResponse 


namespace Surrogate {

    void NonIntrusiveRPCE::Train (
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

        auto NumBasis = 
            BasisFunctions::HermitePolynomials<Z,R,C> (

                NumIndices_, TSet, Dim ()

        );

        auto DenBasis = 
            BasisFunctions::HermitePolynomials<Z,R,C> (

                DenIndices_, TSet, Dim ()

        );

        VectorC NDCoeffs ( NumIndices_.size() + DenIndices_.size() , 0.0 );

        auto nP = NumIndices_.size() / Dim();
        auto nQ = DenIndices_.size() / Dim();

        Eigen::Map<MatrixXC> NumPsi ( 
            NumBasis.data(), NumIndices_.size() / Dim (), nPoints
        );

        // std::cout << NumPsi << std::endl << std::endl;

        Eigen::Map<MatrixXC> DenPsi ( 
            DenBasis.data(), DenIndices_.size() / Dim (), nPoints 
        );

        // std::cout << DenPsi << std::endl << std::endl;

        Eigen::Map<MatrixXC> Coeffs (
            NDCoeffs.data(), 
            ( NumIndices_.size() + DenIndices_.size() ) / Dim (), 
            Dim ()
        );

        Eigen::Map<MatrixXC> Result (
            Response.data(), Dim (), nPoints 
        );

        // std::cout << Result << std::endl;

        NumCoeffs_ = VectorC ( NumIndices_.size(), 0.0 );
        DenCoeffs_ = VectorC ( DenIndices_.size(), 0.0 );

        for ( auto i = 0; i < Dim(); i++ ) {

            DiagMatrixXC M  = Result.row(i).asDiagonal();
            DiagMatrixXC Mh = Result.row(i).conjugate().asDiagonal();

            MatrixXC Aupper ( nP, nP + nQ );
            Aupper << ( NumPsi * NumPsi.transpose() ),
                      ( -NumPsi * M * DenPsi.transpose() );

            // std::cout << Aupper << std::endl << std::endl;

            MatrixXC Alower ( nQ, nP + nQ );
            Alower << ( -DenPsi * Mh * NumPsi.transpose() ), 
                      ( DenPsi * Mh * M * DenPsi.transpose() );

            // std::cout << Alower << std::endl << std::endl;

            MatrixXC A ( nP + nQ, nP + nQ );
            A << Aupper, Alower;

            // std::cout << A << std::endl << std::endl;

            MatrixXC V = A.bdcSvd( Eigen::ComputeFullV ).matrixV();

            VectorXC U = A.bdcSvd( 
                Eigen::ComputeFullV | Eigen::ComputeFullU 
            ).singularValues();

            // if ( i == 1 ) {
            //     // std::cout << V << std::endl; 
            //     // std::cout << std::endl << std::endl;
            //     std::cout << U << std::endl; 
            //
            // }

            VectorXC r = V.col(V.cols() - 1);
            // std::cout << U(V.cols() -1) << std::endl;

            for ( auto j = 0; j < nP; j++ ) {

                NumCoeffs_[i*nP + j] = r(j);

            }

            for ( auto j = 0; j < nQ; j++ ) {
                
                DenCoeffs_[i*nQ+j] = r(nP+j);

            }

        }

    }

} // Non-intrusive RPCE Train 

namespace Surrogate {

    void NonIntrusiveRPCE::PrintCoeffs () const {

        typedef Eigen::Vector<size_t,Eigen::Dynamic> VectorZ;

        Eigen::Map<const VectorXC> NumCoeffs ( 
            NumCoeffs_.data(), NumCoeffs_.size() 
        );

        Eigen::Map<const VectorXC> DenCoeffs ( 
            DenCoeffs_.data(), DenCoeffs_.size() 
        );


        std::cout << NumCoeffs << std::endl;
        std::cout << "-------\n";
        std::cout << DenCoeffs << std::endl;

        for ( auto i = 0; i < NumIndices_.size(); i++ ) {

            std::cout << NumIndices_[i] << std::endl;

        }

    }

} // Non-intrusive PCE PrintCoeffs 

