/**
  * @file GaussSeidel_imp.hpp 
  *
  * @brief implementations of templated Gauss Seidel iterative solver  
  *
  * @author Rezha Adrian Tanuharja
  * Contact: rezha.tanuharja@tum.de / rezhadr@outlook.com 
  */

#ifndef GAUSS_SEIDEL_IMPLEMENTATIONS 
#define GAUSS_SEIDEL_IMPLEMENTATIONS 

#ifndef GAUSS_SEIDEL_DECLARATIONS 
    #include "GaussSeidel.hpp" 
#endif 


namespace linalg {


    template < class Matrix, class Vector, typename T, typename Function >
    GaussSeidel<Matrix,Vector,T,Function>::GaussSeidel ( 
        
        size_t MaxIter, T Tol, const Function& Norm, T omega

    ) : MaxIter_( MaxIter ), Tol_( Tol ), Norm_( Norm ), omega_( omega ) {

        if ( Tol <= 0.0 ) {

            throw std::runtime_error (
                "GaussSeidel Solver: tolerance must be positive"
            );

        }

    }


    template < class Matrix, class Vector, typename T, typename Function >
    void GaussSeidel<Matrix,Vector,T,Function>::UpdateSolution (

        const Matrix& A, const Vector& b, Vector& x 

    ) const {

        if ( A.nCol() != x.size() ||
             A.nRow() != b.size()
        ) {
            throw std::runtime_error (
                "GaussSeidel Solver: sizes don't match"
            );
        }

        Vector e = b - A * x;

        for ( auto i = 0; i < x.size(); i++ ) {
            
            #ifdef LA_COMPLEX
                auto Ax_i = std::complex<T> (0,0);
            #else
                auto Ax_i = 0.0;
            #endif

            for ( auto j = 0; j < x.size(); j++ ) {

                Ax_i += A(i,j) * x[j];

            }

            x[i] += omega_ * ( b[i] - Ax_i ) / A(i,i);

        }

    } // UpdateSolution 


    template < class Matrix, class Vector, typename T, typename Function >
    void GaussSeidel<Matrix,Vector,T,Function>::Solve ( 

        const Matrix& A, const Vector& b, Vector& x

    ) const {

        for ( auto i = 0; i < MaxIter_; i++ ) {

            UpdateSolution ( A, b, x );

            if ( Norm_(b-A*x) < Tol_ ) {
                break;
            }

        }

    } // Solve 


} // linalg 


#endif // GAUSS_SEIDEL_IMPLEMENTATIONS 

